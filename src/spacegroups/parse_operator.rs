use std::iter;

use nom::{IResult, Finish};
use nom::branch::alt;
use nom::character::complete::{char, digit1, one_of, space0};
use nom::combinator::map_opt;
use nom::multi::{many0, separated_list1};
use nom::sequence::{separated_pair, tuple, pair, preceded, delimited};
use num_rational::Ratio;
use num_traits::CheckedAdd;


#[derive(Debug, Eq, Hash, PartialEq)]
enum Axis {
    X,
    Y,
    Z,
    None,
}

type Fraction = Ratio<i32>;
type Term = (Fraction, Axis);
type Coordinate = (Vec<Ratio<i32>>, Ratio<i32>);
type Operator = (Vec<Vec<Ratio<i32>>>, Vec<Ratio<i32>>);


pub fn parse_operator(input: &str) -> Result<(&str, Operator), String>
{
    match delimited(space0, operator, space0)(input).finish() {
        Ok((rest, op)) => Ok((rest, convert_operator(&op))),
        Err(err) => Err(err.to_string()),
    }
}


fn operator(input: &str) -> IResult<&str, Vec<Coordinate>> {
    separated_list1(tuple((space0, char(','), space0)), coordinate)(input)
}


fn coordinate(input: &str) -> IResult<&str, Coordinate> {
    map_opt(
        pair(term, many0(preceded(space0, signed_term))),
        |(first, rest)| collect_terms(first, rest)
    )(input)
}


fn term(input: &str) -> IResult<&str, Term> {
    alt((signed_term, unsigned_term))(input)
}


fn signed_term(input: &str) -> IResult<&str, Term> {
    map_opt(
        separated_pair(sign, space0, unsigned_term),
        |(s, (q, a))| Some((Ratio::from_integer(s) * q, a))
    )(input)
}


fn unsigned_term(input: &str) -> IResult<&str, Term> {
    alt((
        map_opt(
            separated_pair(unsigned_fraction, space0, axis),
            |(q, a)| Some((q, a))
        ),
        map_opt(unsigned_fraction, |q| Some((q, Axis::None))),
        map_opt(axis, |a| Some((Ratio::from_integer(1), a))),
    ))(input)
}


fn unsigned_fraction(input: &str) -> IResult<&str, Fraction> {
    alt((
        map_opt(
            separated_pair(
                integer, tuple((space0, char('/'), space0)), integer
            ),
            |(n, d)| Some(Ratio::new(n as i32, d as i32))
        ),
        map_opt(integer, |n| Some(Ratio::from_integer(n as i32))),
    ))(input)
}


fn sign(input: &str) -> IResult<&str, i32> {
    map_opt(one_of("+-"), map_sign)(input)
}


fn axis(input: &str) -> IResult<&str, Axis> {
    map_opt(one_of("xyz"), map_axis)(input)
}


fn integer(input: &str) -> IResult<&str, u32> {
    map_opt(digit1, map_integer)(input)
}


fn convert_operator(op: &[Coordinate]) -> Operator {
    let mut rows = vec![];
    let mut shift = vec![];

    for (r, s) in op {
        rows.push(r.clone());
        shift.push(*s);
    }

    let n = shift.len();
    let rows = rows.iter()
        .map(|v| v.iter().take(n).copied().collect())
        .collect();

    (rows, shift)
}


fn collect_terms(first: Term, rest: Vec<Term>) -> Option<Coordinate> {
    let mut shift = Ratio::new(0, 1);
    let mut row = vec![Ratio::new(0, 1); 3];

    for (q, a) in iter::once(first).chain(rest) {
        match a {
            Axis::X => row[0] = q.checked_add(&row[0])?,
            Axis::Y => row[1] = q.checked_add(&row[1])?,
            Axis::Z => row[2] = q.checked_add(&row[2])?,
            Axis::None => shift = q.checked_add(&shift)?,
        }
    }

    Some((row, shift))
}


fn map_sign(c: char) -> Option<i32> {
    match c {
        '+' => Some(1),
        '-' => Some(-1),
        _ => None,
    }
}


fn map_axis(c: char) -> Option<Axis> {
    match c {
        'x' => Some(Axis::X),
        'y' => Some(Axis::Y),
        'z' => Some(Axis::Z),
        _ => None,
    }
}


fn map_integer(digits: &str) -> Option<u32> {
    if digits.len() <= 9 {
        Some(digits.chars().fold(0, |n, c| n * 10 + c.to_digit(10).unwrap()))
    } else {
        None
    }
}


#[test]
fn test_parse_sign() {
    assert_eq!(sign("-3x, "), Ok(("3x, ", -1)));
    assert_eq!(sign("+zx, "), Ok(("zx, ", 1)));
    assert!(sign("*, ").is_err());
}


#[test]
fn test_parse_axis() {
    assert_eq!(axis("x, "), Ok((", ", Axis::X)));
    assert_eq!(axis("zx, "), Ok(("x, ", Axis::Z)));
    assert!(axis("a, ").is_err());
}


#[test]
fn test_parse_integer() {
    assert_eq!(integer("123  "), Ok(("  ", 123)));
    assert_eq!(integer("123456789  "), Ok(("  ", 123456789)));
    assert!(integer("x123  ").is_err());
    assert!(integer("1234567890  ").is_err());
}


#[test]
fn test_parse_fraction() {
    assert_eq!(unsigned_fraction("12/3  "), Ok(("  ", Ratio::new(12, 3))));
    assert_eq!(unsigned_fraction("23  "), Ok(("  ", Ratio::new(23, 1))));
    assert!(unsigned_fraction("/3  ").is_err());
}


#[test]
fn test_parse_unsigned_term() {
    assert_eq!(
        unsigned_term("12/3x  "), Ok(("  ", (Ratio::new(12, 3), Axis::X)))
    );
    assert_eq!(
        unsigned_term("12/3  "), Ok(("  ", (Ratio::new(12, 3), Axis::None)))
    );
    assert_eq!(
        unsigned_term("y/  "), Ok(("/  ", (Ratio::new(1, 1), Axis::Y)))
    );
    assert!(unsigned_term("/x  ").is_err());
}


#[test]
fn test_parse_signed_term() {
    assert_eq!(
        signed_term("+12/3x  "), Ok(("  ", (Ratio::new(12, 3), Axis::X)))
    );
    assert_eq!(
        signed_term("-12/3  "), Ok(("  ", (Ratio::new(-12, 3), Axis::None)))
    );
    assert_eq!(
        signed_term("+y/  "), Ok(("/  ", (Ratio::new(1, 1), Axis::Y)))
    );
    assert!(signed_term("4/5x  ").is_err());
    assert!(signed_term("/x  ").is_err());
}


#[test]
fn test_parse_term() {
    assert_eq!(term("12/3x  "), Ok(("  ", (Ratio::new(12, 3), Axis::X))));
    assert_eq!(term("12/3  "), Ok(("  ", (Ratio::new(12, 3), Axis::None))));
    assert_eq!(term("y/  "), Ok(("/  ", (Ratio::new(1, 1), Axis::Y))));
    assert_eq!(term("+12/3 x  "), Ok(("  ", (Ratio::new(12, 3), Axis::X))));
    assert_eq!(
        term("- 12 / 3  "), Ok(("  ", (Ratio::new(-12, 3), Axis::None)))
    );
    assert_eq!(term("+y/  "), Ok(("/  ", (Ratio::new(1, 1), Axis::Y))));
    assert!(term("/x  ").is_err());
}


#[test]
fn test_parse_coordinate() {
    assert_eq!(
        coordinate("x-z"),
        Ok((
            "",
            (
                vec![Ratio::new(1, 1), Ratio::new(0, 1), Ratio::new(-1, 1)],
                Ratio::new(0, 1)
            )
        ))
    );
    assert_eq!(
        coordinate("-1/3 - 2 x + 1 / 2 y - 1/3y + 1/4"),
        Ok((
            "",
            (
                vec![Ratio::new(-2, 1), Ratio::new(1, 6), Ratio::new(0, 1)],
                Ratio::new(-1, 12)
            )
        ))
    );
}


#[test]
fn test_parse_operator_3d() {
    let r = |d| Ratio::new(d, 1);

    assert_eq!(
        parse_operator("  -x, y - x, -z  "),
        Ok((
            "",
            (
                vec![
                    vec![r(-1), r(0), r(0)],
                    vec![r(-1), r(1), r(0)],
                    vec![r(0), r(0), r(-1)],
                ],
                vec![r(0), r(0), r(0)]
            )
        ))
    );
}


#[test]
fn test_parse_operator_2d() {
    let r = |d| Ratio::new(d, 1);

    assert_eq!(
        parse_operator("  1 - x, 2y - x  "),
        Ok((
            "",
            (
                vec![
                    vec![r(-1), r(0)],
                    vec![r(-1), r(2)],
                ],
                vec![r(1), r(0)]
            )
        ))
    );
}
