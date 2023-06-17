use std::io::{BufRead, Read, BufReader};

use super::parse_cgd_line::{parse_cgd_line, Field};


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

    fn add_error(&mut self, msg: &str, lineno: usize) {
        let note = Note::Error(format!("In line {}: {:?}", lineno, msg));
        self.notes.push(note);
    }

    fn add_warning(&mut self, msg: &str, lineno: usize) {
        let note = Note::Warning(format!("In line {}: {:?}", lineno, msg));
        self.notes.push(note);
    }
}


pub fn parse_blocks<T: Read>(input: T) -> Vec<Block> {
    let mut current_block = Block::new();
    let mut blocks = vec![];
    let mut current_key: Option<String> = None;
    let mut in_block = false;

    for (lineno, line) in BufReader::new(input).lines().enumerate() {
        match line {
            Err(e) => current_block.add_error(&format!("{:?}", e), lineno),
            Ok(s) => match parse_cgd_line(&s) {
                Err(e) => current_block.add_error(&e, lineno),
                Ok((rest, fields)) => {
                    if rest.len() > 0 {
                        let msg = format!("unparsed line ending '{:?}'", rest);
                        current_block.add_warning(&msg, lineno);
                    }
                    if fields.len() == 0 {
                        continue
                    }
                }
            }
        }
    }

    blocks
}
