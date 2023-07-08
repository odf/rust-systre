use std::io::{BufRead, Read, BufReader};

use super::parse_cgd_line::{parse_cgd_line, Field};


#[derive(Debug, PartialEq)]
pub enum EntryType {
    Begin, End, Data, Empty, Error
}


#[derive(Debug, PartialEq)]
pub enum Note {
    Error(String),
    Warning(String),
}


#[derive(Debug)]
pub struct Entry {
    entry_type: EntryType,
    key: Option<String>,
    fields: Vec<Field>,
    notes: Vec<Note>,
    line_number: usize,
    content: Option<String>,
}


impl Entry {
    fn new(
        entry_type: EntryType,
        line_number: usize,
        content: Option<String>,
        key: Option<String>,
        fields: Vec<Field>
    )
        -> Self
    {
        Entry { entry_type, line_number, content, key, fields, notes: vec![] }
    }

    fn add_error(&mut self, msg: &str) {
        self.notes.push(Note::Error(msg.to_string()));
    }

    fn add_warning(&mut self, msg: &str) {
        self.notes.push(Note::Warning(msg.to_string()));
    }
}


#[derive(Debug)]
pub struct Block {
    block_type: String,
    lineno_start: usize,
    lineno_end: usize,
    entries: Vec<Entry>,
}


impl Block {
    fn new() -> Self {
        Block {
            block_type: String::new(),
            lineno_start: 0,
            lineno_end: 0,
            entries: vec![],
        }
    }
}


pub fn parse_blocks<T: Read>(input: T) -> Vec<Block> {
    let mut blocks = vec![];
    let mut current_block = Block::new();
    let mut current_key: Option<String> = None;
    let mut in_block = false;

    for mut entry in parse_lines(input) {
        match entry.entry_type {
            EntryType::Empty => {
                if entry.notes.len() > 0 {
                    current_block.entries.push(entry);
                }
            },
            EntryType::Error => {
                current_block.entries.push(entry);
            },
            EntryType::Begin | EntryType::End => {
                // should not happen
            },
            EntryType::Data => {
                if let Some(s) = entry.key {
                    let new_key = s.to_lowercase();
                    entry.key = Some(new_key.clone());

                    if new_key == "end" {
                        entry.entry_type = EntryType::End;
                        if !in_block {
                            current_block.lineno_start = entry.line_number;
                            current_block.block_type ="--void--".to_string();
                            entry.add_error("block has no type or content");
                        }
                        current_block.lineno_end = entry.line_number;
                        if entry.fields.len() > 1 {
                            entry.add_warning("text after 'end' ignored");
                        }
                        current_block.entries.push(entry);
                        blocks.push(current_block);
                        current_block = Block::new();
                        in_block = false;
                    } else if !in_block {
                        entry.entry_type = EntryType::Begin;
                        current_block.lineno_start = entry.line_number;
                        current_block.block_type = new_key;
                        in_block = true;
                        current_key = None;
                        if entry.fields.len() > 1 {
                            entry.add_warning("text after block type ignored");
                        }
                        current_block.entries.push(entry);
                    } else {
                        current_key = Some(new_key.clone());
                        current_block.entries.push(entry);
                    }
                } else {
                    entry.key = current_key.clone();
                    if !in_block {
                        entry.add_error("data found before block start");
                    }
                    current_block.entries.push(entry);
                }
            },
        }
    }

    if in_block {
        let last_entry = current_block.entries.last_mut().unwrap();

        current_block.lineno_end = last_entry.line_number;
        last_entry.add_error("final block is missing an 'end' statement");
        blocks.push(current_block);
    }

    blocks
}


pub fn parse_lines<T: Read>(input: T) -> Vec<Entry> {
    let mut all_entries = vec![];
    let mut lineno = 0;

    for line in BufReader::new(input).lines() {
        let mut current_entry = Entry::new(
            EntryType::Data, lineno, None, None, vec![]
        );

        match line {
            Err(e) => {
                current_entry.entry_type = EntryType::Error;
                current_entry.add_error(&format!("{:?}", e));
            },
            Ok(s) => match parse_cgd_line(&s) {
                Err(e) => {
                    current_entry.entry_type = EntryType::Error;
                    current_entry.content = Some(s);
                    current_entry.add_error(&e);
                },
                Ok((rest, fields)) => {
                    current_entry.content = Some(s.clone());
                    if rest.len() > 0 {
                        let msg = format!("unparsed line ending '{}'", rest);
                        current_entry.add_error(&msg);
                    }

                    if fields.len() == 0 {
                        current_entry.entry_type = EntryType::Empty;
                    } else if let Field::Key(s) = &fields[0] {
                        current_entry.key = Some(s.to_lowercase());
                        current_entry.fields = fields[1..].to_vec();
                    } else {
                        current_entry.fields = fields;
                    }
                }
            }
        }

        all_entries.push(current_entry);
        lineno += 1;
    }

    all_entries
}


#[test]
fn test_parse_cgd_block() {
    let input =
r#"
FIRST
    Data 1 two 27/25
         3 four +27.25E-3
    Name "Gavrog"
    Data "some more" "s"a
END

SECOND
    456
END

123

THIRD
"#;

    let blocks = parse_blocks(input.as_bytes());

    assert_eq!(blocks.len(), 3);

    let first = &blocks[0];
    eprintln!("{:?}", first);
    assert_eq!(first.block_type, "first");
    assert_eq!(first.lineno_start, 1);
    assert_eq!(first.lineno_end, 6);

    assert_eq!(first.entries.len(), 6);
    assert_eq!(first.entries[0].entry_type, EntryType::Begin);
    assert_eq!(first.entries[1].entry_type, EntryType::Data);
    assert_eq!(first.entries[1].key, Some(String::from("data")));
    assert_eq!(first.entries[2].key, Some(String::from("data")));
    assert_eq!(first.entries[3].key, Some(String::from("name")));
    assert_eq!(first.entries[4].key, Some(String::from("data")));
    assert_eq!(first.entries[5].entry_type, EntryType::End);

    assert_eq!(first.entries[4].notes.len(), 1);
    assert_eq!(
        first.entries[4].notes[0],
        Note::Error("unparsed line ending 'a'".to_string())
    );

    let second = &blocks[1];
    assert_eq!(second.block_type, "second");
    assert_eq!(second.lineno_start, 8);
    assert_eq!(second.lineno_end, 10);

    assert_eq!(second.entries.len(), 3);
    assert_eq!(second.entries[1].notes.len(), 0);
    assert_eq!(second.entries[1].key, None);

    let third = &blocks[2];
    assert_eq!(third.block_type, "third");
    assert_eq!(third.lineno_start, 14);
    assert_eq!(third.lineno_end, 14);

    assert_eq!(third.entries.len(), 2);
    assert_eq!(third.entries[0].notes.len(), 1);
    assert_eq!(
        third.entries[0].notes[0],
        Note::Error("data found before block start".to_string())
    );
    assert_eq!(third.entries[1].notes.len(), 1);
    assert_eq!(
        third.entries[1].notes[0],
        Note::Error("final block is missing an 'end' statement".to_string())
    );
}
