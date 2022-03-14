pub fn quick_sort<T: PartialOrd>(v: &mut Vec<T>, p: usize, r: usize) {
    if p < r {
        let q = partition(v, p, r);
        quick_sort(v, p, q);
        quick_sort(v, q + 1, r);
    }
}

fn partition<T: PartialOrd>(v: &mut Vec<T>, p: usize, r: usize) -> usize {
    let mut i = p as i32 - 1;
    let mut j = r as i32 + 1;

    loop {
        while {
            j -= 1;
            v[j as usize] > v[p]
        } {}
        while {
            i += 1;
            v[i as usize] < v[p]
        } {}
        if i < j {
            // rustだと入れ替えをちょっと頑張る。unsafeでやってもいいけど
            // let p1: *mut T = &mut v[i as usize];
            // let p2: *mut T = &mut v[j as usize];
            // unsafe {
            //     p1.swap(p2);
            // }
            let (v1, v2) = v.split_at_mut(j as usize); // [0,j)と[j, v.len())に分ける
            std::mem::swap(&mut v1[i as usize], &mut v2[0]);
        } else {
            return j as usize;
        }
    }
}
