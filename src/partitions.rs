use std::cell::UnsafeCell;
use std::hash::Hash;
use std::collections::HashMap;


struct PartitionImpl<T> {
    rank: HashMap<T, usize>,
    parent: HashMap<T, T>,
}

impl<T> PartitionImpl<T> where T: Copy + Eq + Hash {
    fn new() -> Self {
        PartitionImpl { rank: HashMap::new(), parent: HashMap::new() }
    }

    fn find(&mut self, x: &T) -> T {
        let mut root = x;
        while self.parent.contains_key(root) {
            root = self.parent.get(root).unwrap();
        }
        let root = *root;

        let mut x = *x;
        while x != root {
            let t = x;
            x = *self.parent.get(&x).unwrap();
            self.parent.insert(t, root);
        }

        root
    }

    fn unite(&mut self, x: &T, y: &T) {
        let x0 = self.find(x);
        let y0 = self.find(y);

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


impl<T> Partition<T> where T: Copy + Eq + Hash {
    pub fn new() -> Self {
        Partition { _impl: UnsafeCell::new(PartitionImpl::new())}
    }

    pub fn find(&self, x: &T) -> T {
        unsafe { (*self._impl.get()).find(x) }
    }

    pub fn unite(&mut self, x: &T, y: &T) {
        unsafe { (*self._impl.get()).unite(x, y) };
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
}
