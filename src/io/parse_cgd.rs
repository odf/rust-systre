use nom::bytes::complete::{take_till1, take_until1};
use nom::character::complete::{char, digit1, space0};
use nom::character::is_space;
use nom::number::complete::double;
use nom::sequence::{separated_pair, tuple, preceded, delimited};
use nom::combinator::map_opt;
use nom::branch::alt;
use nom::IResult;
use num_rational::Ratio;


enum Field {
    Fraction(Ratio<i32>),
    Double(f64),
    Int(i32),
    Name(String),
}


fn field(input: &str) -> IResult<&str, Field> {
    alt((
        map_opt(fraction, |f| Some(Field::Fraction(f))),
        map_opt(double, |f| Some(Field::Double(f))),
        map_opt(integer, |f| Some(Field::Int(f))),
        map_opt(quoted_string, |f| Some(Field::Name(f.to_string()))),
        map_opt(name, |f| Some(Field::Name(f.to_string()))),
    ))(input)
}


fn fraction(input: &str) -> IResult<&str, Ratio<i32>> {
    alt((
        map_opt(
            separated_pair(
                integer,
                tuple((space0, char('/'), space0)),
                unsigned_integer
            ),
            |(n, d)| Some(Ratio::new(n, d as i32))
        ),
        map_opt(unsigned_integer, |n| Some(Ratio::from_integer(n as i32))),
    ))(input)
}


fn integer(input: &str) -> IResult<&str, i32> {
    alt((
        map_opt(unsigned_integer, |n| Some(n as i32)),
        map_opt(preceded(char('-'), unsigned_integer), |n| Some(-(n as i32)))
    ))(input)
}


fn quoted_string(input: &str) -> IResult<&str, &str> {
    delimited(char('"'), take_until1("\""), char('"'))(input)
}


fn name(input: &str) -> IResult<&str, &str> {
    take_till1(|c| c == '"' || is_space(c as u8))(input)
}


fn unsigned_integer(input: &str) -> IResult<&str, u32> {
    map_opt(digit1, map_integer)(input)
}


fn map_integer(digits: &str) -> Option<u32> {
    if digits.len() <= 9 {
        Some(digits.chars().fold(0, |n, c| n * 10 + c.to_digit(10).unwrap()))
    } else {
        None
    }
}


#[test]
fn test_parse_integer() {
    assert_eq!(integer("-123  "), Ok(("  ", -123)));
    assert_eq!(integer("123456789  "), Ok(("  ", 123456789)));
    assert!(integer("x123  ").is_err());
    assert!(integer("1234567890  ").is_err());
}


#[test]
fn test_parse_fraction() {
    assert_eq!(fraction("12 / 3  "), Ok(("  ", Ratio::new(12, 3))));
    assert_eq!(fraction("-12/23  "), Ok(("  ", Ratio::new(-12, 23))));
    assert_eq!(fraction("23  "), Ok(("  ", Ratio::new(23, 1))));
    assert!(fraction("/3  ").is_err());
}


#[test]
fn test_parse_quoted_string() {
    assert_eq!(quoted_string("\" \"  "), Ok(("  ", " ")));
    assert_eq!(quoted_string("\"abc \"  "), Ok(("  ", "abc ")));
    assert!(quoted_string("\"abc").is_err());
}
