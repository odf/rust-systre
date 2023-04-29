use nom::IResult;
use nom::branch::alt;
use nom::character::complete::{char, digit1, one_of, space0};
use nom::combinator::map_opt;
use nom::multi::{many0, separated_list1};
use nom::sequence::{separated_pair, tuple, pair, preceded};


#[derive(Debug, PartialEq)]
enum Axis {
    X,
    Y,
    Z,
    None,
}

type Fraction = (i32, u32);
type Term = (Fraction, Axis);
type Coordinate = Vec<Term>;
type Operator = Vec<Coordinate>;


fn operator(input: &str) -> IResult<&str, Operator> {
    separated_list1(tuple((space0, char(','), space0)), coordinate)(input)
}


fn coordinate(input: &str) -> IResult<&str, Coordinate> {
    map_opt(
        pair(term, many0(preceded(space0, signed_term))),
        |(first, rest)| Some(collect_terms(first, rest))
    )(input)
}


fn term(input: &str) -> IResult<&str, Term> {
    alt((signed_term, unsigned_term))(input)
}


fn signed_term(input: &str) -> IResult<&str, Term> {
    map_opt(
        separated_pair(sign, space0, unsigned_term),
        |(s, ((n, d), a))| Some(((s * n, d), a))
    )(input)
}


fn unsigned_term(input: &str) -> IResult<&str, Term> {
    alt((
        map_opt(
            separated_pair(unsigned_fraction, space0, axis),
            |(q, a)| Some((q, a))
        ),
        map_opt(unsigned_fraction, |q| Some((q, Axis::None))),
        map_opt(axis, |a| Some(((1, 1), a))),
    ))(input)
}


fn unsigned_fraction(input: &str) -> IResult<&str, Fraction> {
    alt((
        map_opt(
            separated_pair(
                integer, tuple((space0, char('/'), space0)), integer
            ),
            |(n, d)| Some((n as i32, d))
        ),
        map_opt(integer, |n| Some((n as i32, 1))),
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


fn collect_terms(first: Term, rest: Vec<Term>) -> Coordinate {
    let mut terms = vec![first];
    for t in rest {
        terms.push(t);
    }
    terms
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
    Some(digits.chars().fold(0, |n, c| n * 10 + c.to_digit(10).unwrap()))
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
    assert!(integer("x123  ").is_err());
}


#[test]
fn test_parse_fraction() {
    assert_eq!(unsigned_fraction("12/3  "), Ok(("  ", (12, 3))));
    assert_eq!(unsigned_fraction("23  "), Ok(("  ", (23, 1))));
    assert!(unsigned_fraction("/3  ").is_err());
}


#[test]
fn test_parse_unsigned_term() {
    assert_eq!(unsigned_term("12/3x  "), Ok(("  ", ((12, 3), Axis::X))));
    assert_eq!(unsigned_term("12/3  "), Ok(("  ", ((12, 3), Axis::None))));
    assert_eq!(unsigned_term("y/  "), Ok(("/  ", ((1, 1), Axis::Y))));
    assert!(unsigned_term("/x  ").is_err());
}


#[test]
fn test_parse_signed_term() {
    assert_eq!(signed_term("+12/3x  "), Ok(("  ", ((12, 3), Axis::X))));
    assert_eq!(signed_term("-12/3  "), Ok(("  ", ((-12, 3), Axis::None))));
    assert_eq!(signed_term("+y/  "), Ok(("/  ", ((1, 1), Axis::Y))));
    assert!(signed_term("4/5x  ").is_err());
    assert!(signed_term("/x  ").is_err());
}


#[test]
fn test_parse_term() {
    assert_eq!(term("12/3x  "), Ok(("  ", ((12, 3), Axis::X))));
    assert_eq!(term("12/3  "), Ok(("  ", ((12, 3), Axis::None))));
    assert_eq!(term("y/  "), Ok(("/  ", ((1, 1), Axis::Y))));
    assert_eq!(term("+12/3 x  "), Ok(("  ", ((12, 3), Axis::X))));
    assert_eq!(term("- 12 / 3  "), Ok(("  ", ((-12, 3), Axis::None))));
    assert_eq!(term("+y/  "), Ok(("/  ", ((1, 1), Axis::Y))));
    assert!(term("/x  ").is_err());
}


#[test]
fn test_parse_coordinate() {
    assert_eq!(
        coordinate("x-z"),
        Ok(("", vec![((1, 1), Axis::X), ((-1, 1), Axis::Z)]))
    );
    assert_eq!(
        coordinate("-1/3 - 2 x + 1 / 2 y"),
        Ok((
            "",
            vec![((-1, 3), Axis::None), ((-2, 1), Axis::X), ((1, 2), Axis::Y)]
        ))
    );
}


#[test]
fn test_parse_operator() {
    assert_eq!(
        operator("-x, y - x, -z"),
        Ok((
            "",
            vec![
                vec![((-1, 1), Axis::X)],
                vec![((1, 1), Axis::Y), ((-1, 1), Axis::X)],
                vec![((-1, 1), Axis::Z)],
            ]
        ))
    );
}
