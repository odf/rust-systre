use nom::IResult;
use nom::branch::alt;
use nom::character::complete::{char, digit1, one_of};
use nom::combinator::map_opt;
use nom::sequence::separated_pair;


fn map_axis(c: char) -> Option<u8> {
    match c {
        'x' => Some(0),
        'y' => Some(1),
        'z' => Some(2),
        _ => None,
    }
}


fn map_integer(digits: &str) -> Option<u32> {
    Some(digits.chars().fold(0, |n, c| n * 10 + c.to_digit(10).unwrap()))
}


fn axis(input: &str) -> IResult<&str, u8> {
    map_opt(one_of("xyz"), map_axis)(input)
}


fn integer(input: &str) -> IResult<&str, u32> {
    map_opt(digit1, map_integer)(input)
}


fn fraction(input: &str) -> IResult<&str, (u32, u32)> {
    let quotient = separated_pair(integer, char('/'), integer);
    let simple = map_opt(integer, |n| Some((n, 1)));

    alt((quotient, simple))(input)
}


#[test]
fn test_parser_axis() {
    assert_eq!(axis("x, "), Ok((", ", 0)));
    assert_eq!(axis("zx, "), Ok(("x, ", 2)));
    assert!(axis("a, ").is_err());
}


#[test]
fn test_parser_integer() {
    assert_eq!(integer("123  "), Ok(("  ", 123)));
    assert!(integer("x123  ").is_err());
}


#[test]
fn test_parser_fraction() {
    assert_eq!(fraction("12/3  "), Ok(("  ", (12, 3))));
    assert_eq!(fraction("23  "), Ok(("  ", (23, 1))));
    assert!(fraction("/3  ").is_err());
}
