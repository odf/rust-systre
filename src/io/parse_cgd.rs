use nom::character::complete::{char, digit1, space0};
use nom::sequence::{separated_pair, tuple, preceded};
use nom::combinator::map_opt;
use nom::branch::alt;
use nom::IResult;
use num_rational::Ratio;


type Fraction = Ratio<i32>;


fn fraction(input: &str) -> IResult<&str, Fraction> {
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
        integer,
        map_opt(preceded(char('-'), integer), |n| Some(-n))
    ))(input)
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
