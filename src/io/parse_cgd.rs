use nom::{IResult, Finish};
use nom::bytes::complete::{take_till1, take_until1, tag};
use nom::character::complete::{char, digit1, space0, space1};
use nom::character::is_space;
use nom::multi::separated_list1;
use nom::number::complete::double;
use nom::sequence::{separated_pair, preceded, delimited, pair};
use nom::combinator::{map_opt, rest};
use nom::branch::alt;
use num_rational::Ratio;


#[derive(Debug, PartialEq)]
pub enum Field {
    Fraction(Ratio<i32>),
    Int(i32),
    Double(f64),
    Name(String),
}


pub fn parse_cgd_line(input: &str) -> Result<(&str, Vec<Field>), String> {
    match line(input).finish() {
        Ok((rest, fields)) => Ok((rest, fields)),
        Err(err) => Err(err.to_string()),
    }
}


fn line(input: &str) -> IResult<&str, Vec<Field>> {
    alt((
        map_opt(
            separated_pair(separated_list1(space1, field), space0, comment),
            |(fields, _)| Some(fields)
        ),
        separated_list1(space1, field),
    ))(input)
}


fn comment(input: &str) -> IResult<&str, ()> {
    map_opt(pair(tag("#"), rest), |_| Some(()))(input)
}


fn field(input: &str) -> IResult<&str, Field> {
    alt((
        map_opt(fraction, |f| Some(Field::Fraction(f))),
        map_opt(
            double,
            |f| if f == (f as i32) as f64 {
                Some(Field::Int(f as i32))
            } else {
                Some(Field::Double(f))
            }
        ),
        map_opt(quoted_string, |f| Some(Field::Name(f.to_string()))),
        map_opt(name, |f| Some(Field::Name(f.to_string()))),
    ))(input)
}


fn fraction(input: &str) -> IResult<&str, Ratio<i32>> {
    map_opt(
        separated_pair(integer, char('/'), unsigned_integer),
        |(n, d)| Some(Ratio::new(n, d as i32))
    )(input)
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
    take_till1(|c| c == '"' || c == '#' || is_space(c as u8))(input)
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
fn test_parse_cgd_integer() {
    assert_eq!(integer("-123  "), Ok(("  ", -123)));
    assert_eq!(integer("123456789  "), Ok(("  ", 123456789)));
    assert!(integer("x123  ").is_err());
    assert!(integer("1234567890  ").is_err());
}


#[test]
fn test_parse_cgd_fraction() {
    assert_eq!(fraction("12/3  "), Ok(("  ", Ratio::new(12, 3))));
    assert_eq!(fraction("-12/23  "), Ok(("  ", Ratio::new(-12, 23))));
    assert!(fraction("/3  ").is_err());
}


#[test]
fn test_parse_cgd_quoted_string() {
    assert_eq!(quoted_string("\" \"  "), Ok(("  ", " ")));
    assert_eq!(quoted_string("\"abc \"  "), Ok(("  ", "abc ")));
    assert!(quoted_string("\"abc").is_err());
}


#[test]
fn test_parse_cgd_field() {
    assert_eq!(field("-123  "), Ok(("  ", Field::Int(-123))));
    assert_eq!(field("-123.4  "), Ok(("  ", Field::Double(-123.4))));
    assert_eq!(
        field("12/34  "),
        Ok(("  ", Field::Fraction(Ratio::new(12, 34))))
    );
    assert_eq!(
        field("\"a s d\"  "),
        Ok(("  ", Field::Name("a s d".to_string())))
    );
    assert_eq!(field("xyz  "), Ok(("  ", Field::Name("xyz".to_string()))));
}


#[test]
fn test_parse_cgd_line() {
    assert_eq!(
        parse_cgd_line("ATOM \"Si1\" 4 .5 -0.3 2.1 # This is a comment "),
        Ok(("",
            vec![
                Field::Name("ATOM".to_string()),
                Field::Name("Si1".to_string()),
                Field::Int(4),
                Field::Double(0.5),
                Field::Double(-0.3),
                Field::Double(2.1),
            ]
        ))
    );
    assert_eq!(
        parse_cgd_line("asdf \""),
        Ok((" \"", vec![Field::Name("asdf".to_string())]))
    );
}
