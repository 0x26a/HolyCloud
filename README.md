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

let sphere = vec![[5.000000000000001, 5.000000000000001, 7.0710678118654755], [7.0710678118654755, 7.0710678118654755, 6.123233995736766e-16], [5.000000000000001, 5.000000000000001, -7.071067811865475], [4.329780281177467e-16, 7.0710678118654755, 7.0710678118654755], [6.123233995736766e-16, 10.0, 6.123233995736766e-16], [4.329780281177467e-16, 7.0710678118654755, -7.071067811865475], [-5.0, 5.000000000000001, 7.0710678118654755], [-7.071067811865475, 7.0710678118654755, 6.123233995736766e-16], [-5.0, 5.000000000000001, -7.071067811865475], [-7.0710678118654755, 8.659560562354934e-16, 7.0710678118654755], [-10.0, 1.2246467991473533e-15, 6.123233995736766e-16], [-7.0710678118654755, 8.659560562354934e-16, -7.071067811865475], [-5.000000000000002, -5.0, 7.0710678118654755], [-7.071067811865477, -7.071067811865475, 6.123233995736766e-16], [-5.000000000000002, -5.0, -7.071067811865475], [-1.2989340843532398e-15, -7.0710678118654755, 7.0710678118654755], [-1.8369701987210296e-15, -10.0, 6.123233995736766e-16], [-1.2989340843532398e-15, -7.0710678118654755, -7.071067811865475], [4.999999999999999, -5.000000000000002, 7.0710678118654755], [7.071067811865474, -7.071067811865477, 6.123233995736766e-16], [4.999999999999999, -5.000000000000002, -7.071067811865475]];

// VietorisRips::<d, n> : compute the nth first homology modules in d-dimensional space
let mut vr = VietorisRips::<3, 3>::init(sphere);

// vr.analyse(ring, radius limit, radius step)
let data = vr.analyze(Ring::Z, 22.0, 0.5);
for modules in data{
  println!("\nHOMOLOGICAL DATA ( persistence={} )", modules.end - module.start);
  for i in 0..3{
    println!("\tdim H_{} = {} | T(H_{}) = {:?}", i, modules.ranks[i],i, modules.torsions[i]);
  }
}
```

## Output
```
HHOMOLOGICAL DATA ( persistence=7.5 )
        dim H_0 = 9 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []
        dim H_2 = 0 | T(H_2) = []

HOMOLOGICAL DATA ( persistence=2.5 )
        dim H_0 = 7 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []
        dim H_2 = 0 | T(H_2) = []

HOMOLOGICAL DATA ( persistence=0.5 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 6 | T(H_1) = []
        dim H_2 = 0 | T(H_2) = []

HOMOLOGICAL DATA ( persistence=1 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 20 | T(H_1) = []
        dim H_2 = 0 | T(H_2) = []

HOMOLOGICAL DATA ( persistence=3 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 3 | T(H_1) = []
        dim H_2 = 0 | T(H_2) = []

HOMOLOGICAL DATA ( persistence=3 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 1 | T(H_1) = []
        dim H_2 = 0 | T(H_2) = []

HOMOLOGICAL DATA ( persistence=4 )
        dim H_0 = 1 | T(H_0) = []
        dim H_1 = 0 | T(H_1) = []
        dim H_2 = 1 | T(H_2) = []
```
We see that the two most persistent topological features are:
* The homology of the discrete point cloud
    - $$H_0(V_\epsilon)=\mathbb{Z}^9$$
    - $$H_1(V_\epsilon)=0$$
    - $$H_2(V_\epsilon)=0$$
* The homology of the sphere $$S^2$$
  - $$H_0(V_\epsilon)=\mathbb{Z}$$
  - $$H_1(V_\epsilon)=0$$
  - $$H_2(V_\epsilon)=\mathbb{Z}$$
 
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
