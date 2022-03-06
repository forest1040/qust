use srand::{Rand, RngSource};

pub struct Xor128 {
    x: u64,
    y: u64,
    z: u64,
    w: u64,
}

impl Xor128 {
    pub fn from_seed(seed: i64) -> Xor128 {
        let mut r: Rand<_> = Rand::new(RngSource::new(seed));
        let ignore_first = 40;
        let mut res = Xor128 {
            x: r.uint64(),
            y: r.uint64(),
            z: r.uint64(),
            w: r.uint64(),
        };
        for _ in 0..ignore_first {
            res.next();
        }
        res
    }

    pub fn next(&mut self) -> u64 {
        let t = self.x ^ (self.x << 11);
        self.x = self.y;
        self.y = self.z;
        self.z = self.w;
        self.w = (self.w ^ (self.w >> 19)) ^ (t ^ (t >> 8));
        self.w
    }

    pub fn random_uniform(&mut self) -> f64 {
        self.next() as f64 / std::u64::MAX as f64
    }

    pub fn random_normal(&mut self) -> f64 {
        // return sqrt(-1.0 * log(random_uniform(state))) *
        //     sin(2.0 * M_PI * random_uniform(state));
        (-1.0 * self.random_uniform().log(std::f64::consts::E)).sqrt()
            * (2.0 * std::f64::consts::PI * self.random_uniform()).sin()
    }
}
