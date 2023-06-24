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
    entries: Vec<(String, Vec<Field>)>,
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

    fn add_entry(&mut self, key: &str, fields: &[Field]) {
        self.entries.push((key.to_string(), fields.iter().cloned().collect()));
    }
}


pub fn parse_blocks<T: Read>(input: T) -> Vec<Block> {
    let mut current_block = Block::new();
    let mut blocks = vec![];
    let mut current_key: Option<String> = None;
    let mut in_block = false;
    let mut lineno = 0;

    for line in BufReader::new(input).lines() {
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
                        // nothing to do for this line
                    } else if let Field::Key(s) = &fields[0] {
                        let new_key = s.to_lowercase();

                        if new_key == "end" {
                            if !in_block {
                                current_block.lineno_start = lineno;
                                current_block.block_type =
                                    "--invalid--".to_string();
                                let msg = "block has no type or content";
                                current_block.add_error(msg, lineno, None);
                            }
                            current_block.lineno_end = lineno;
                            if fields.len() > 1 {
                                let msg = "text after 'end' keyword ignored";
                                current_block.add_warning(
                                    msg, lineno, Some(&s)
                                );
                            }
                            blocks.push(current_block);

                            current_block = Block::new();
                            in_block = false;
                        } else if !in_block {
                            current_block.lineno_start = lineno;
                            current_block.block_type = new_key;
                            in_block = true;
                            if fields.len() > 1 {
                                let msg = "text after block type ignored";
                                current_block.add_warning(
                                    msg, lineno, Some(&s)
                                );
                            }
                        } else {
                            current_key = Some(new_key.clone());
                            current_block.add_entry(&new_key, &fields[1..]);
                        }
                    } else if let Some(key) = &current_key {
                        current_block.add_entry(&key, &fields);
                    } else {
                        let msg = "block data is not preceded by a keyword";
                        current_block.add_error(msg, lineno, Some(&s));
                    }
                }
            }
        }

        lineno += 1;
    }

    if in_block {
        let msg = "final block is missing and 'end' statement";
        current_block.lineno_end = lineno;
        current_block.add_error(msg, 0, None);
        blocks.push(current_block);
    }

    blocks
}
