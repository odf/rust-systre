use std::io::stdin;

use rust_systre::spacegroups::parse_table::parse_space_group_table;


fn main() {
    parse_space_group_table(stdin());
}
