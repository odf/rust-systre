use std::io::{BufRead, Read, BufReader};

use crate::io::parse_cgd_line::Field;


pub enum Note {
    Error(String),
    Warning(String),
}


pub struct Block {
    block_type: String,
    lineno_start: u32,
    lineno_end: u32,
    entries: Vec<Vec<Field>>,
    notes: Vec<Note>,
}


impl Block {
    fn new() -> Self {
        Block {
            block_type: String::new(),
            lineno_start: 0,
            lineno_end: 0,
            entries: vec![],
            notes: vec![]
        }
    }
}


pub fn parse_blocks<T: Read>(input: T) -> Vec<Block> {
    let mut blocks = vec![];
    let mut current_block = Block::new();
    let mut in_block = false;

    for (lineno, line) in BufReader::new(input).lines().enumerate() {

    }

    blocks
}
