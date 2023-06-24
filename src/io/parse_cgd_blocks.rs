use std::io::{BufRead, Read, BufReader};

use super::parse_cgd_line::{parse_cgd_line, Field};


pub enum Note {
    Error(String),
    Warning(String),
}


pub struct Context {
    line_number: usize,
    line_content: Option<String>,
}


pub struct Block {
    block_type: String,
    lineno_start: usize,
    lineno_end: usize,
    entries: Vec<Vec<Field>>,
    notes: Vec<(Note, Context)>,
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

    fn add_note(
        &mut self, note: Note, line_number: usize, line: Option<&str>
    ) {
        let line_content = line.and_then(|s| Some(s.to_string()));
        self.notes.push((note, Context { line_number, line_content }));
    }

    fn add_error(
        &mut self, msg: &str, line_number: usize, line: Option<&str>
    ) {
        self.add_note(Note::Error(msg.to_string()), line_number, line);
    }

    fn add_warning(
        &mut self, msg: &str, line_number: usize, line: Option<&str>
    ) {
        self.add_note(Note::Warning(msg.to_string()), line_number, line);
    }
}


pub fn parse_blocks<T: Read>(input: T) -> Vec<Block> {
    let mut current_block = Block::new();
    let mut blocks = vec![];
    let mut current_key: Option<String> = None;
    let mut in_block = false;

    for (lineno, line) in BufReader::new(input).lines().enumerate() {
        match line {
            Err(e) => {
                current_block.add_error(&format!("{:?}", e), lineno, None)
            },
            Ok(s) => match parse_cgd_line(&s) {
                Err(e) => current_block.add_error(&e, lineno, Some(&s)),
                Ok((rest, fields)) => {
                    if rest.len() > 0 {
                        let msg = format!("unparsed line ending '{:?}'", rest);
                        current_block.add_warning(&msg, lineno, Some(&s));
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
