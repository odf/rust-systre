use std::cell::UnsafeCell;
use std::hash::Hash;
use std::collections::HashMap;


struct PartitionImpl<T> {
    elements: Vec<T>,
    index: HashMap<T, usize>,
    rank: HashMap<usize, usize>,
    parent: HashMap<usize, usize>,
}

impl<T> PartitionImpl<T> where T: Clone + Eq + Hash {
    fn new() -> Self {
        PartitionImpl {
            elements: vec![],
            index: HashMap::new(),
            rank: HashMap::new(),
            parent: HashMap::new(),
        }
    }

    fn get_index(&mut self, a: &T) -> usize {
        if self.index.contains_key(&a) {
            *self.index.get(&a).unwrap()
        } else {
            let i = self.elements.len();
            self.elements.push(a.clone());
            self.index.insert(a.clone(), i);
            i
        }
    }

    fn find(&mut self, a: &T) -> T {
        let mut x = self.get_index(a);
        let mut root = x;

        while self.parent.contains_key(&root) {
            root = *self.parent.get(&root).unwrap();
        }

        while x != root {
            let t = x;
            x = *self.parent.get(&x).unwrap();
            self.parent.insert(t, root);
        }

        self.elements[root].clone()
    }

    fn unite(&mut self, a: &T, b: &T) {
        let a0 = self.find(a);
        let b0 = self.find(b);

        let x0 = *self.index.get(&a0).unwrap();
        let y0 = *self.index.get(&b0).unwrap();

        if x0 != y0 {
            let rx = *self.rank.get(&x0).unwrap_or(&0);
            let ry = *self.rank.get(&y0).unwrap_or(&0);

            if rx < ry {
                self.parent.insert(x0, y0);
            } else {
                if rx == ry {
                    self.rank.insert(x0, rx + 1);
                }
                self.parent.insert(y0, x0);
            }
        }
    }
}


pub struct Partition<T> {
    _impl: UnsafeCell<PartitionImpl<T>>,
}


impl<T> Partition<T> where T: Clone + Eq + Hash {
    pub fn new() -> Self {
        Partition { _impl: UnsafeCell::new(PartitionImpl::new())}
    }

    pub fn find(&self, x: &T) -> T {
        unsafe { (*self._impl.get()).find(x) }
    }

    pub fn unite(&mut self, x: &T, y: &T) {
        unsafe { (*self._impl.get()).unite(x, y) };
    }

    pub fn classes(&self, elms: &[T]) -> Vec<Vec<T>> {
        let mut class_for_rep = HashMap::new();
        let mut classes = vec![];

        for e in elms {
            let rep = self.find(e);
            if let Some(cl) = class_for_rep.get(&rep) {
                let class: &mut Vec<_> = &mut classes[*cl];
                class.push(e.clone());
            } else {
                class_for_rep.insert(rep, classes.len());
                classes.push(vec![e.clone()]);
            }
        }

        classes
    }
}


#[test]
pub fn test_partition() {
    let p = {
        let mut p = Partition::new();
        for (a, b) in [(1, 2), (3, 4), (5, 6), (7, 8), (2, 3), (1, 6)] {
            p.unite(&a, &b);
        }
        p
    };

    let test = HashMap::from([
        (0, 0),
        (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1),
        (7, 2), (8, 2),
        (9, 3),
    ]);

    for a in 0..=9 {
        for b in 0..=9 {
            if test[&a] == test[&b] {
                assert_eq!(p.find(&a), p.find(&b));
            } else {
                assert_ne!(p.find(&a), p.find(&b));
            }
        }
    }

    let elms: Vec<_> = (0..=9).collect();
    let cl = p.classes(&elms);
    assert_eq!(cl, vec![vec![0], vec![1, 2, 3, 4, 5, 6], vec![7, 8], vec![9]]);
}
