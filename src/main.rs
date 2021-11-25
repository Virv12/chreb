use num_rational::Rational32;

#[derive(Debug)]
struct Item<'a> {
    s: &'a str,
    n: u32,
    mol: bool,
    arr: bool,
}

fn build_linsys(v: &[Item]) -> (usize, Vec<Rational32>) {
    let n = v.iter().filter(|&i| i.mol).count();
    let mut m = Vec::new();
    let mut s = Vec::new();
    let mut k = 0;

    let mut mul = 1;

    for i in v {
        let p = match s.iter().position(|&x| x == i.s) {
            Some(p) => p,
            None => {
                s.push(i.s);
                m.resize(m.len() + n, Rational32::new(0, 1));
                s.len() - 1
            }
        };

        k += i.mol as usize;
        mul -= 2 * i.arr as i32;
        m[p * n + k - 1] += i.n as i32 * mul;
    }

    debug_assert_eq!(k, n);

    (n, m)
}

fn solve_linsys(n: usize, m: &mut Vec<Rational32>) -> Vec<usize> {
    assert!(m.len() % n == 0);
    let k = m.len() / n;
    let mut a = Vec::new();
    let mut j = 0;

    for i in 0..n {
        if let Some(e) = (j..k)
            .filter(|x| m[x * n + i] != Rational32::new(0, 1))
            .next()
        {
            for t in 0..n {
                m.swap(j * n + t, e * n + t);
            }
        } else {
            a.push(i);
            continue;
        }

        for e in 0..k {
            let f = if e == j {
                Rational32::new(1, 1) - m[j * n + i].recip()
            } else {
                m[e * n + i] / m[j * n + i]
            };
            for t in 0..n {
                let q = f * m[j * n + t];
                m[e * n + t] -= q;
            }
        }

        j += 1;
    }

    m.resize(j * n, Rational32::new(0, 1));

    a
}

fn parse(mut s: &str) -> Option<Vec<Item>> {
    fn parse_name(s: &str) -> Option<(&str, &str)> {
        let mut it = s.chars();
        if !it.next()?.is_ascii_uppercase() {
            return None;
        }

        let i = 1 + it.take_while(char::is_ascii_lowercase).count();
        Some((&s[i..], &s[..i]))
    }

    fn parse_int(s: &str) -> Option<(&str, u32)> {
        let it = s.chars();
        let i = it.take_while(char::is_ascii_digit).count();
        if i == 0 {
            return None;
        }
        Some((&s[i..], s[..i].parse().ok()?))
    }

    let mut v = Vec::new();

    let mut arr_total = false;
    let mut arr = false;
    loop {
        let mut mol = true;
        while let Some((s2, name)) = parse_name(s) {
            let (s2, n) = parse_int(s2).unwrap_or((s2, 1));
            v.push(Item {
                s: name,
                n,
                mol,
                arr,
            });
            s = s2;
            mol = false;
            arr = false;
        }

        let i = s.chars().take_while(char::is_ascii_whitespace).count();
        s = &s[i..];

        let mut it = s.chars();
        s = match it.next() {
            Some('+') => &s[1..],
            Some('-') => {
                if it.next() != Some('>') {
                    return None;
                }

                if arr_total {
                    return None;
                }

                arr = true;
                arr_total = true;
                &s[2..]
            }
            None => break,
            _ => return None,
        };

        let i = s.chars().take_while(char::is_ascii_whitespace).count();
        s = &s[i..];
    }

    if !arr_total {
        return None;
    }

    Some(v)
}

fn main() {
    let mut buf = String::new();
    let stdin = std::io::stdin();
    stdin.read_line(&mut buf).unwrap();

    let v = parse(&buf).expect("Input is malformed");
    let (n, mut ls) = build_linsys(&v);
    let dof = solve_linsys(n, &mut ls);

    for &i in &dof {
        let p = (0..ls.len() / n)
            .map(|x| *ls[x * n + i].denom())
            .fold(1, num_integer::lcm);

        let mut k = 0;
        for m in 0..n {
            let c = if m == i {
                p
            } else if dof.contains(&m) {
                0
            } else {
                let e = (0..ls.len() / n)
                    .filter(|&x| ls[x * n + m] == Rational32::new(1, 1))
                    .next();
                (-ls[e.unwrap() * n + i] * p).to_integer()
            };

            if k != 0 {
                if v[k].arr {
                    print!(" -> ");
                } else {
                    print!(" + ");
                }
            }
            if c != 1 {
                print!("{}", c);
            }
            loop {
                print!("{}", v[k].s);
                if v[k].n != 1 {
                    print!("{}", v[k].n);
                }
                k += 1;
                if k == v.len() || v[k].mol {
                    break;
                }
            }
        }
        print!("\n");
    }
}
