# HolyCloud
A Rust library to compute persistent homology with torsion of point clouds in 2D and 3D.

## Warnings
This library was written for learning purpose and may be unstable.


## Details
* K-d trees are used to build the Vietoris-Rips complex.
* torsion of homology $$\mathbb{Z}$$-modules is computed with a naive implementation of Smith Normal Form.
  - [Modules over PID, part II. Smith Normal Form.](https://m-ershov.github.io/7752_Spring2010/lecture8.pdf)
  - [Calculating Homology of a Simplicial Complex Using Smith Normal Form](https://eric-bunch.github.io/blog/calculating_homology_of_simplicial_complex)
* SVD is computed with the nalgebra Rust binding for Lapack (FORTRAN).

## Example

Persistent homology of a sphere-like point cloud in 3D:
```rust
use holycloud::topology::Ring;
use holycloud::persistence::VietorisRips;

let two_circles = vec![[4.504844339512096, 2.1694186955877908, 10.0], [3.1174490092936677, 3.909157412340149, 10.0], [1.1126046697815721, 4.874639560909118, 10.0], [-1.1126046697815717, 4.874639560909118, 10.0], [-3.1174490092936673, 3.9091574123401496, 10.0], [-4.504844339512095, 2.169418695587791, 10.0], [-5.0, 6.123233995736766e-16, 10.0], [-4.504844339512096, -2.16941869558779, 10.0], [-3.1174490092936686, -3.9091574123401482, 10.0], [-1.112604669781573, -4.874639560909118, 10.0], [1.1126046697815712, -4.874639560909118, 10.0], [3.117449009293667, -3.9091574123401496, 10.0], [4.504844339512095, -2.1694186955877917, 10.0], [4.504844339512096, 2.1694186955877908, 0.0], [3.1174490092936677, 3.909157412340149, 0.0], [1.1126046697815721, 4.874639560909118, 0.0], [-1.1126046697815717, 4.874639560909118, 0.0], [-3.1174490092936673, 3.9091574123401496, 0.0], [-4.504844339512095, 2.169418695587791, 0.0], [-5.0, 6.123233995736766e-16, 0.0], [-4.504844339512096, -2.16941869558779, 0.0], [-3.1174490092936686, -3.9091574123401482, 0.0], [-1.112604669781573, -4.874639560909118, 0.0], [1.1126046697815712, -4.874639560909118, 0.0], [3.117449009293667, -3.9091574123401496, 0.0], [4.504844339512095, -2.1694186955877917, 0.0]];

// VietorisRips::<d, n> : compute the nth first homology modules in d-dimensional space
let mut vr = VietorisRips::<3, 2>::init(two_circles);

// vr.analyse(ring, radius limit, radius step)
let data = vr.analyze(Ring::R, 20.0, 0.5);
for modules in data{
  println!("\nHOMOLOGICAL DATA ( persistence={} )", modules.end - module.start);
  for i in 0..2{
    println!("\tdim H_{} = {} | T(H_{}) = {:?}", i, modules.ranks[i],i, modules.torsions[i]);
  }
}
```

## Output
```
HOMOLOGICAL DATA ( persistence=2.5 )
        dim H_0 = 26 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 22 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 10 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=1 )
        dim H_0 = 2 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=5.5 )
        dim H_0 = 2 | T(H_0) = []
        dim H_1 = 2 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 2 | T(H_0) = []
        dim H_1 = 4 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 14 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=1.5 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 12 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 10 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 4 | T(H_1) = []

HOMOLOGICAL DATA ( persistence=7 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []
```
We see that the three most persistent topological features are:
* The homology of two circles
    - $$H_0(V_\epsilon)=\mathbb{Z}^2$$
    - $$H_1(V_\epsilon)=\mathbb{Z}^2$$
* The homology of the initial discrete point clouds
  - $$H_0(V_\epsilon)=\mathbb{Z}^{26}$$
  - $$H_1(V_\epsilon)=0$$
* The homology of the ending simplicial complex (homotopic to a point)
  - $$H_0(V_\epsilon)=\mathbb{Z}$$
  - $$H_1(V_\epsilon)=0$$

 
## More Examples
You can compute the $$\mathbb{Z}$$ or $$\mathbb{R}$$ simplicial homology of a few samples.
```rust
use holycloud::topology::SimplicialComplex;
use holycloud::topology::samples::{S1,T2,RP2};

let rp2 = SimplicialComplex::const_from(&RP2, Ring::Z);

// hx.0 = rank, hx.1 = torsions
let (h0, h1, h2) = rp2.H(0), rp2.H(1), rp2.H(2);

assert_eq!(h0.0, 1) // Z
assert_eq!(h0.1, vec![]) // 0

assert_eq!(h1.0, 0) // 0
assert_eq!(h1.1, vec![2]) // Z/2Z

assert_eq!(h2.0, 0) // 0
assert_eq!(h2.1, vec![]) // 0

```

And of your own simplicial complexes.
```rust
// from a constant array of simplicies
const X: [&const_Simplex;3] = [
	&const_Simplex{value: &[1], dimension: 0},
	&const_Simplex{value: &[2], dimension: 0},
	&const_Simplex{value: &[1,2], dimension: 1}, // every simplicies must be oriented in increasing order.
];
let x = SimplicialComplex::const_from(&X, Ring::R);
```
```rust
// from a vector of simplicies
let X: Vec<Simplex> = vec![
  Simplex::from(vec![1]),
  Simplex::from(vec![2]),
  Simplex::from(vec![1,2]),
];
let x = SimplicialComplex::from(X, Ring::R);
```

Or maybe just compute their Betti numbers.
```rust
// 1st Betti number of the torus
let x = SimplicialComplex::const_from(&T2, Ring::R);
assert_eq!(x.b(1), 2);
```
