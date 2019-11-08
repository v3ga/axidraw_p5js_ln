

//  
// https://github.com/fogleman/ln/blob/master/ln/util.go
//
const Radians = (degrees) => degrees * Math.PI / 180;
const Degrees = (radians) => radians * 180 / Math.PI;
const Median = (items) => {
	let n = items.length;
	if (n == 0) {
		return 0;
	} else if (n%2 == 1) {
		return items[n/2];
	} else {
		const a = items[n/2-1];
		const b = items[n/2];
		return (a + b) / 2;
	}
}
const INF = 1e9;
const EPS = 1e-9;
//
// https://github.com/fogleman/ln/blob/master/ln/axis.go
//
const AxisNone = -1;
const AxisX = 0;
const AxisY = 1;
const AxisZ = 2;
//  
// https://github.com/fogleman/ln/blob/master/ln/vector.go 
//
class Vector {
    constructor(x, y, z) {
    	this.X = x; this.Y = y; this.Z = z;
    }
    RandomUnitVector() {
        let x, y, z;
    	do {
    		x = Math.ranomd()*2 - 1;
    		y = Math.ranomd()*2 - 1;
    		z = Math.ranomd()*2 - 1;
    	} while (x*x+y*y+z*z > 1);
    	return new Vector(x, y, z).Normalize();
    }
    Length() { return Math.sqrt(this.X*this.X + this.Y*this.Y + this.Z*this.Z); }
    Distance(b) { return this.Sub(b).Length(); }
    LengthSquared() { return this.X*this.X + this.Y*this.Y + this.Z*this.Z; }
    DistanceSquared(b){ return this.Sub(b).LengthSquared(); }
    Dot(b) { return this.X*b.X + this.Y*b.Y + this.Z*b.Z; }
    Cross(b) { 
        return new Vector(  this.Y*b.Z - this.Z*b.Y, 
                            this.Z*b.X - this.X*b.Z, 
                            this.X*b.Y - this.Y*b.X); }
    Normalize() { return this.MulScalar(1/this.Length()); }
    Add(b) { return new Vector(this.X + b.X, this.Y + b.Y, this.Z + b.Z); }
    Sub(b) { return new Vector(this.X - b.X, this.Y - b.Y, this.Z - b.Z); }
    Mul(b) { return new Vector(this.X * b.X, this.Y * b.Y, this.Z * b.Z); }
    Div(b) { return new Vector(this.X / b.X, this.Y / b.Y, this.Z / b.Z); }
    AddScalar(b) { return new Vector(this.X + b, this.Y + b, this.Z + b); }
    SubScalar(b) { return new Vector(this.X - b, this.Y - b, this.Z - b); }
    MulScalar(b) { return new Vector(this.X * b, this.Y * b, this.Z * b); }
    DivScalar(b) { return new Vector(this.X / b, this.Y / b, this.Z / b); }
    Min(b) { return new Vector(Math.min(this.X, b.X), Math.min(this.Y, b.Y), Math.min(this.Z, b.Z)); }
    Max(b) { return new Vector(Math.max(this.X, b.X), Math.max(this.Y, b.Y), Math.max(this.Z, b.Z)); }
    MinAxis() {
    	const x = Math.abs(this.X), y = Math.abs(this.Y), z = Math.abs(this.Z);
    	if (x <= y && x <= z) {
    		return new Vector(1, 0, 0);
    	} else if (y <= x && y <= z) {
    		return new Vector(0, 1, 0);
    	}
    	return new Vector(0, 0, 1);
    }
    MinComponent() { return Math.min(Math.min(this.X, this.Y), this.Z); }
    SegmentDistance(v, w) {
    	const l2 = v.DistanceSquared(w)
    	if (l2 == 0) {
    		return this.Distance(v);
    	}
    	const t = this.Sub(v).Dot(w.Sub(v)) / l2;
    	if (t < 0) {
    		return this.Distance(v);
    	}
    	if (t > 1){
    		return this.Distance(w);
    	}
    	return v.Add(w.Sub(v).MulScalar(t)).Distance(this);
    }
}

//
// https://github.com/fogleman/ln/blob/master/ln/matrix.go
//
class Matrix {
    constructor(x00, x01, x02, x03,
            	x10, x11, x12, x13,
	            x20, x21, x22, x23,
	            x30, x31, x32, x33) {
	               	this.x00=x00, this.x01=x01, this.x02=x02, this.x03=x03;
	               	this.x10=x10, this.x11=x11, this.x12=x12, this.x13=x13;
	               	this.x20=x20, this.x21=x21, this.x22=x22, this.x23=x23;
	               	this.x30=x30, this.x31=x31, this.x32=x32, this.x33=x33;
	            }
    Translate(v) {
    	return Translate(v).Mul(this);
    }
    Scale(v) {
    	return Scale(v).Mul(this);
    }
    Rotate(v, a) {
    	return Rotate(v, a).Mul(this);
    }
    Frustum(l, r, b, t, n, f) {
    	return Frustum(l, r, b, t, n, f).Mul(this);
    }
    Orthographic(l, r, b, t, n, f) {
    	return Orthographic(l, r, b, t, n, f).Mul(this);
    }
    Perspective(fovy, aspect, near, far) {
    	return Perspective(fovy, aspect, near, far).Mul(this);
    }
    Mul(b) {
    	const m = new Matrix();
    	m.x00 = this.x00*b.x00 + this.x01*b.x10 + this.x02*b.x20 + this.x03*b.x30;
    	m.x10 = this.x10*b.x00 + this.x11*b.x10 + this.x12*b.x20 + this.x13*b.x30;
    	m.x20 = this.x20*b.x00 + this.x21*b.x10 + this.x22*b.x20 + this.x23*b.x30;
    	m.x30 = this.x30*b.x00 + this.x31*b.x10 + this.x32*b.x20 + this.x33*b.x30;
    	m.x01 = this.x00*b.x01 + this.x01*b.x11 + this.x02*b.x21 + this.x03*b.x31;
    	m.x11 = this.x10*b.x01 + this.x11*b.x11 + this.x12*b.x21 + this.x13*b.x31;
    	m.x21 = this.x20*b.x01 + this.x21*b.x11 + this.x22*b.x21 + this.x23*b.x31;
    	m.x31 = this.x30*b.x01 + this.x31*b.x11 + this.x32*b.x21 + this.x33*b.x31;
    	m.x02 = this.x00*b.x02 + this.x01*b.x12 + this.x02*b.x22 + this.x03*b.x32;
    	m.x12 = this.x10*b.x02 + this.x11*b.x12 + this.x12*b.x22 + this.x13*b.x32;
    	m.x22 = this.x20*b.x02 + this.x21*b.x12 + this.x22*b.x22 + this.x23*b.x32;
    	m.x32 = this.x30*b.x02 + this.x31*b.x12 + this.x32*b.x22 + this.x33*b.x32;
    	m.x03 = this.x00*b.x03 + this.x01*b.x13 + this.x02*b.x23 + this.x03*b.x33;
    	m.x13 = this.x10*b.x03 + this.x11*b.x13 + this.x12*b.x23 + this.x13*b.x33;
    	m.x23 = this.x20*b.x03 + this.x21*b.x13 + this.x22*b.x23 + this.x23*b.x33;
    	m.x33 = this.x30*b.x03 + this.x31*b.x13 + this.x32*b.x23 + this.x33*b.x33;
    	return m;
    }
    MulPosition(b) {
    	const x = this.x00*b.X + this.x01*b.Y + this.x02*b.Z + this.x03;
    	const y = this.x10*b.X + this.x11*b.Y + this.x12*b.Z + this.x13;
    	const z = this.x20*b.X + this.x21*b.Y + this.x22*b.Z + this.x23;
    	return new Vector(x, y, z);
    }
    MulPositionW(b) {
    	const x = this.x00*b.X + this.x01*b.Y + this.x02*b.Z + this.x03;
    	const y = this.x10*b.X + this.x11*b.Y + this.x12*b.Z + this.x13;
    	const z = this.x20*b.X + this.x21*b.Y + this.x22*b.Z + this.x23;
    	const w = this.x30*b.X + this.x31*b.Y + this.x32*b.Z + this.x33;
    	return new Vector(x / w, y / w, z / w);
    }
    MulDirection(b) {
    	const x = this.x00*b.X + this.x01*b.Y + this.x02*b.Z;
    	const y = this.x10*b.X + this.x11*b.Y + this.x12*b.Z;
    	const z = this.x20*b.X + this.x21*b.Y + this.x22*b.Z;
    	return new Vector(x, y, z).Normalize();
    }
    MulRay(b) {
    	return new Ray(this.MulPosition(b.Origin), this.MulDirection(b.Direction));
    }
    MulBox(box) {
    	// http://dev.theomader.com/transform-bounding-boxes/
    	const r = new Vector(this.x00, this.x10, this.x20);
    	const u = new Vector(this.x01, this.x11, this.x21);
    	const b = new Vector(this.x02, this.x12, this.x22);
    	const t = new Vector(this.x03, this.x13, this.x23);
    	let xa = r.MulScalar(box.Min.X);
    	let xb = r.MulScalar(box.Max.X);
    	let ya = u.MulScalar(box.Min.Y);
    	let yb = u.MulScalar(box.Max.Y);
    	let za = b.MulScalar(box.Min.Z);
    	let zb = b.MulScalar(box.Max.Z);
    	const min = xa.Min(xb).Add(ya.Min(yb)).Add(za.Min(zb)).Add(t);
    	const max = xa.Max(xb).Add(ya.Max(yb)).Add(za.Max(zb)).Add(t);
    	return new Box(min, max);
    }
    Transpose() {
    	return new Matrix(
    		this.x00, this.x10, this.x20, this.x30,
    		this.x01, this.x11, this.x21, this.x31,
    		this.x02, this.x12, this.x22, this.x32,
    		this.x03, this.x13, this.x23, this.x33);
    }
    Determinant() {
    	return (this.x00*this.x11*this.x22*this.x33 - this.x00*this.x11*this.x23*this.x32 +
    		this.x00*this.x12*this.x23*this.x31 - this.x00*this.x12*this.x21*this.x33 +
    		this.x00*this.x13*this.x21*this.x32 - this.x00*this.x13*this.x22*this.x31 -
    		this.x01*this.x12*this.x23*this.x30 + this.x01*this.x12*this.x20*this.x33 -
    		this.x01*this.x13*this.x20*this.x32 + this.x01*this.x13*this.x22*this.x30 -
    		this.x01*this.x10*this.x22*this.x33 + this.x01*this.x10*this.x23*this.x32 +
    		this.x02*this.x13*this.x20*this.x31 - this.x02*this.x13*this.x21*this.x30 +
    		this.x02*this.x10*this.x21*this.x33 - this.x02*this.x10*this.x23*this.x31 +
    		this.x02*this.x11*this.x23*this.x30 - this.x02*this.x11*this.x20*this.x33 -
    		this.x03*this.x10*this.x21*this.x32 + this.x03*this.x10*this.x22*this.x31 -
    		this.x03*this.x11*this.x22*this.x30 + this.x03*this.x11*this.x20*this.x32 -
    		this.x03*this.x12*this.x20*this.x31 + this.x03*this.x12*this.x21*this.x30);
    }
    Inverse() {
    	const m = new Matrix();
    	const d = this.Determinant();
    	m.x00 = (this.x12*this.x23*this.x31 - this.x13*this.x22*this.x31 + this.x13*this.x21*this.x32 - this.x11*this.x23*this.x32 - this.x12*this.x21*this.x33 + this.x11*this.x22*this.x33) / d;
    	m.x01 = (this.x03*this.x22*this.x31 - this.x02*this.x23*this.x31 - this.x03*this.x21*this.x32 + this.x01*this.x23*this.x32 + this.x02*this.x21*this.x33 - this.x01*this.x22*this.x33) / d;
    	m.x02 = (this.x02*this.x13*this.x31 - this.x03*this.x12*this.x31 + this.x03*this.x11*this.x32 - this.x01*this.x13*this.x32 - this.x02*this.x11*this.x33 + this.x01*this.x12*this.x33) / d;
    	m.x03 = (this.x03*this.x12*this.x21 - this.x02*this.x13*this.x21 - this.x03*this.x11*this.x22 + this.x01*this.x13*this.x22 + this.x02*this.x11*this.x23 - this.x01*this.x12*this.x23) / d;
    	m.x10 = (this.x13*this.x22*this.x30 - this.x12*this.x23*this.x30 - this.x13*this.x20*this.x32 + this.x10*this.x23*this.x32 + this.x12*this.x20*this.x33 - this.x10*this.x22*this.x33) / d;
    	m.x11 = (this.x02*this.x23*this.x30 - this.x03*this.x22*this.x30 + this.x03*this.x20*this.x32 - this.x00*this.x23*this.x32 - this.x02*this.x20*this.x33 + this.x00*this.x22*this.x33) / d;
    	m.x12 = (this.x03*this.x12*this.x30 - this.x02*this.x13*this.x30 - this.x03*this.x10*this.x32 + this.x00*this.x13*this.x32 + this.x02*this.x10*this.x33 - this.x00*this.x12*this.x33) / d;
    	m.x13 = (this.x02*this.x13*this.x20 - this.x03*this.x12*this.x20 + this.x03*this.x10*this.x22 - this.x00*this.x13*this.x22 - this.x02*this.x10*this.x23 + this.x00*this.x12*this.x23) / d;
    	m.x20 = (this.x11*this.x23*this.x30 - this.x13*this.x21*this.x30 + this.x13*this.x20*this.x31 - this.x10*this.x23*this.x31 - this.x11*this.x20*this.x33 + this.x10*this.x21*this.x33) / d;
    	m.x21 = (this.x03*this.x21*this.x30 - this.x01*this.x23*this.x30 - this.x03*this.x20*this.x31 + this.x00*this.x23*this.x31 + this.x01*this.x20*this.x33 - this.x00*this.x21*this.x33) / d;
    	m.x22 = (this.x01*this.x13*this.x30 - this.x03*this.x11*this.x30 + this.x03*this.x10*this.x31 - this.x00*this.x13*this.x31 - this.x01*this.x10*this.x33 + this.x00*this.x11*this.x33) / d;
    	m.x23 = (this.x03*this.x11*this.x20 - this.x01*this.x13*this.x20 - this.x03*this.x10*this.x21 + this.x00*this.x13*this.x21 + this.x01*this.x10*this.x23 - this.x00*this.x11*this.x23) / d;
    	m.x30 = (this.x12*this.x21*this.x30 - this.x11*this.x22*this.x30 - this.x12*this.x20*this.x31 + this.x10*this.x22*this.x31 + this.x11*this.x20*this.x32 - this.x10*this.x21*this.x32) / d;
    	m.x31 = (this.x01*this.x22*this.x30 - this.x02*this.x21*this.x30 + this.x02*this.x20*this.x31 - this.x00*this.x22*this.x31 - this.x01*this.x20*this.x32 + this.x00*this.x21*this.x32) / d;
    	m.x32 = (this.x02*this.x11*this.x30 - this.x01*this.x12*this.x30 - this.x02*this.x10*this.x31 + this.x00*this.x12*this.x31 + this.x01*this.x10*this.x32 - this.x00*this.x11*this.x32) / d;
    	m.x33 = (this.x01*this.x12*this.x20 - this.x02*this.x11*this.x20 + this.x02*this.x10*this.x21 - this.x00*this.x12*this.x21 - this.x01*this.x10*this.x22 + this.x00*this.x11*this.x22) / d;
    	return m;
    }
}
function Identity() {
	return new Matrix(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);
}
function Translate(v) {
	return new Matrix(
		1, 0, 0, v.X,
		0, 1, 0, v.Y,
		0, 0, 1, v.Z,
		0, 0, 0, 1);
}
function Scale(v) {
	return new Matrix(
		v.X, 0, 0, 0,
		0, v.Y, 0, 0,
		0, 0, v.Z, 0,
		0, 0, 0, 1);
}
function Rotate(v, a) {
	v = v.Normalize();
	const s = Math.sin(a);
	const c = Math.cos(a);
	const m = 1 - c;
	return new Matrix(
		m*v.X*v.X + c, m*v.X*v.Y + v.Z*s, m*v.Z*v.X - v.Y*s, 0,
		m*v.X*v.Y - v.Z*s, m*v.Y*v.Y + c, m*v.Y*v.Z + v.X*s, 0,
		m*v.Z*v.X + v.Y*s, m*v.Y*v.Z - v.X*s, m*v.Z*v.Z + c, 0,
		0, 0, 0, 1);
}
function Frustum(l, r, b, t, n, f) {
	const t1 = 2 * n;
	const t2 = r - l;
	const t3 = t - b;
	const t4 = f - n;
	return new Matrix(
		t1 / t2, 0, (r + l) / t2, 0,
		0, t1 / t3, (t + b) / t3, 0,
		0, 0, (-f - n) / t4, (-t1 * f) / t4,
		0, 0, -1, 0);
}
function Orthographic(l, r, b, t, n, f) {
	return new Matrix(
		2 / (r - l), 0, 0, -(r + l) / (r - l),
		0, 2 / (t - b), 0, -(t + b) / (t - b),
		0, 0, -2 / (f - n), -(f + n) / (f - n),
		0, 0, 0, 1);
}
function Perspective(fovy, aspect, near, far) {
	const ymax = near * Math.tan(fovy*Math.PI/360);
	const xmax = ymax * aspect;
	return Frustum(-xmax, xmax, -ymax, ymax, near, far);
}
function LookAt(eye, center, up) {
	up = up.Normalize();
	const f = center.Sub(eye).Normalize();
	const s = f.Cross(up).Normalize();
	const u = s.Cross(f).Normalize();
	const m = new Matrix(
		s.X, u.X, -f.X, eye.X,
		s.Y, u.Y, -f.Y, eye.Y,
		s.Z, u.Z, -f.Z, eye.Z,
		0, 0, 0, 1,
	);
	return m.Inverse();
}
//  
// https://github.com/fogleman/ln/blob/master/ln/tree.go
//
class Tree {
    constructor(shapes) {
    	this.Box = BoxForShapes(shapes);
    	this.Root = new Node(shapes);
    	this.Root.Split(0);
    }
    Intersect(r) {
    	const i = this.Box.Intersect(r);
    	if (i.t2 < i.t1 || i.t2 <= 0) {
    		return NoHit;
    	}
    	return this.Root.Intersect(r, i.t1, i.t2);
    }
}

class Node {
    constructor(shapes) {
    	this.Axis =  AxisNone
    	this.Point = 0;
    	this.Shapes = shapes;
    	this.Left  = false;
    	this.Right = false;
    }
    Intersect(r, tmin, tmax) {
    	let tsplit = 0;
    	let leftFirst  = false;
    	switch (this.Axis) {
        	case AxisNone:
        		return this.IntersectShapes(r);
        	case AxisX:
        		tsplit = (this.Point - r.Origin.X) / r.Direction.X;
        		leftFirst = (r.Origin.X < this.Point) || (r.Origin.X == this.Point && r.Direction.X <= 0);
        	case AxisY:
        		tsplit = (this.Point - r.Origin.Y) / r.Direction.Y;
        		leftFirst = (r.Origin.Y < this.Point) || (r.Origin.Y == this.Point && r.Direction.Y <= 0);
        	case AxisZ:
        		tsplit = (this.Point - r.Origin.Z) / r.Direction.Z;
        		leftFirst = (r.Origin.Z < this.Point) || (r.Origin.Z == this.Point && r.Direction.Z <= 0);
    	}
    	let first, second;
    	if (leftFirst) {
    		first = this.Left;
    		second = this.Right;
    	} else {
    		first = this.Right;
    		second = this.Left;
    	}
    	if (tsplit > tmax || tsplit <= 0) {
    		return first.Intersect(r, tmin, tmax);
    	} else if (tsplit < tmin) {
    		return second.Intersect(r, tmin, tmax);
    	} else {
    		const h1 = first.Intersect(r, tmin, tsplit);
    		if (h1.T <= tsplit) {
    			return h1;
    		}
    		const h2 = second.Intersect(r, tsplit, Math.min(tmax, h1.T));
    		if (h1.T <= h2.T) {
    			return h1;
    		} else {
    			return h2;
    		}
    	}
    }
    IntersectShapes(r) {
    	let hit = NoHit;
    	for (const _ in this.Shapes) {
    	    const shape = this.Shapes[_];
    		const h = shape.Intersect(r);
    		if (h.T < hit.T ){
    			hit = h;
    		}
    	}
    	return hit;
    }
    PartitionScore(axis, point) {
    	let left = 0, right = 0;
    	for (const _ in this.Shapes) {
    	    const shape = this.Shapes[_];
    		const box = shape.BoundingBox();
    		const p = box.Partition(axis, point);
    		if (p.left) {
    			left++;
    		}
    		if (p.right) {
    			right++;
    		}
    	}
    	if (left >= right) {
    		return left;
    	} else {
    		return right;
    	}
    }
    Partition(size, axis, point) {
    	const left = [], right = [];
    	for (const _ in this.Shapes) {
    	    const shape = this.Shapes[_];
    		const box = shape.BoundingBox();
    		const p = box.Partition(axis, point);
    		if (p.left) {
    			left.push(shape);
    		}
    		if (p.right) {
    			right.push(shape);
    		}
    	}
    	return { left, right };
    }
    Split(depth) {
    	if (this.Shapes.length < 8) {
    		return;
    	}
    	const xs = [], ys = [], zs = [];
    	for (const _ in this.Shapes) {
    	    const shape = this.Shapes[_];
    		const box = shape.BoundingBox();
    		xs.push(box.Min.X);
    		xs.push(box.Max.X);
    		ys.push(box.Min.Y);
    		ys.push(box.Max.Y);
    		zs.push(box.Min.Z);
    		zs.push(box.Max.Z);
    	}
    	xs.sort();
    	ys.sort();
    	zs.sort();
    	const mx = Median(xs), my = Median(ys), mz = Median(zs);
    	let best = Math.floor(this.Shapes.length * 0.85);
    	let bestAxis = AxisNone;
    	let bestPoint = 0.0;
    	const sx = this.PartitionScore(AxisX, mx);
    	if (sx < best) {
    		best = sx;
    		bestAxis = AxisX;
    		bestPoint = mx;
    	}
    	const sy = this.PartitionScore(AxisY, my);
    	if (sy < best) {
    		best = sy;
    		bestAxis = AxisY;
    		bestPoint = my;
    	}
    	const sz = this.PartitionScore(AxisZ, mz);
    	if (sz < best) {
    		best = sz;
    		bestAxis = AxisZ;
    		bestPoint = mz;
    	}
    	if (bestAxis == AxisNone ){
    		return;
    	}
    	const p = this.Partition(best, bestAxis, bestPoint)
    	this.Axis = bestAxis;
    	this.Point = bestPoint;
    	this.Left = new Node(p.left);
    	this.Right = new Node(p.right);
    	this.Left.Split(depth + 1);
    	this.Right.Split(depth + 1);
    	this.Shapes = false; // only needed at leaf nodes
    }
}
//  
// https://github.com/fogleman/ln/blob/master/ln/shape.go
//
class TransformedShape{
    constructor(shape, matrix) {
        this.Shape = shape;
        this.Matrix = matrix;
        this.Inverse = matrix.Inverse();
    }
    Compile() {
        this.Shape.Compile();
    }
    BoundingBox() {
    	return this.Matrix.MulBox(this.Shape.BoundingBox());
    }
    Contains(v, f) {
    	return this.Shape.Contains(this.Inverse.MulPosition(v), f);
    }
    Intersect(r ) {
    	return this.Shape.Intersect(this.Inverse.MulRay(r));
    }
    Paths() {
    	return this.Shape.Paths().Transform(this.Matrix);
    }
}
//  
// https://github.com/fogleman/ln/blob/master/ln/scene.go
//
class Scene {
    constructor() {
        this.Shapes = [];
        this.Tree = false;
    }
    Compile() {
    	for (const _ in this.Shapes) {
    	    const shape = this.Shapes[_];
    		shape.Compile()
    	}
    	if (!this.Tree) {
    		this.Tree = new Tree(this.Shapes);
    	}
    }
    Add(shapes) {
    	if(Array.isArray(shapes)) {
    	    this.Shapes = [...this.Shapes, ...shapes];
    	} else {
    	    this.Shapes.push(shapes);
    	}
    }
    Intersect(r) {
    	return this.Tree.Intersect(r);
    }
    Visible(eye, point) {
    	const v = eye.Sub(point);
    	const r = new Ray(point, v.Normalize());
    	const hit = this.Intersect(r);
    	return (hit.T >= v.Length());
    }
    Paths(p = -1) {
        if (p >= 0) {
            if (p < this.Shapes.length) {
                return this.Shapes[p].Paths();
            } else {
                return false;
            }
        } else {
        	let result = new Paths();
        	for (const _ in this.Shapes) {
        	    const shape = this.Shapes[_];
        		result.Append(shape.Paths().v);
        	}
        	return result;
        }
    }
    Render(eye, center, up, width, height, fovy, near, far, step, p = -1) {
    	const aspect = width / height;
    	let matrix = LookAt(eye, center, up);
    	matrix = matrix.Perspective(fovy, aspect, near, far);
    	return this.RenderWithMatrix(matrix, eye, width, height, step, p);
    }
    RenderWithMatrix(matrix, eye, width, height, step, p = -1) {
        if (p <= 0) {
        	this.Compile();
        }
    	let paths = this.Paths(p);
    	if (!paths) {
    	    return false;
    	}
    	if (step > 0) {
    		paths = paths.Chop(step);
    	}
    	paths = paths.Filter(new ClipFilter(matrix, eye, this));
    	if (step > 0) {
    		paths = paths.Simplify(1e-6);
    	}
    	matrix = Scale( new Vector(-width / 2, -height / 2, 0));
    	paths = paths.Transform(matrix);
    	return paths;
    }
}
//
// https://github.com/fogleman/ln/blob/master/ln/ray.go
//
class Ray {
    constructor(origin, direction) {
	    this.Origin = origin;
	    this.Direction = direction;
    }
    Position(t) {
	    return this.Origin.Add(this.Direction.MulScalar(t));
    }
}
//
// https://github.com/fogleman/ln/blob/master/ln/path.go
//
class Path {
    constructor(v = []) {
        this.v = v;
    }
    Append(v) {
        if (Array.isArray(v)) {
            this.v = this.v.concat(v);
        } else {
            this.v.push(v);
        }
    }
    BoundingBox() {
    	let box = new Box(this.v[0], this.v[0]);
    	for (const _ in this.v) {
    	    const v = this.v[_];
    		box = box.Extend(new Box(v, v));
    	}
    	return box
    }
    Transform(matrix) {
    	const result = new Path();
    	for (const _ in this.v) {
    	    const v = this.v[_];
    		result.Append(matrix.MulPosition(v));
    	}
    	return result;
    }
    Chop(step) {
    	const result = new Path();
    	for (let i = 0; i < this.v.length-1; i++) {
    		const a = this.v[i];
    		const b = this.v[i+1];
    		const v = b.Sub(a);
    		const l = v.Length();
    		if (i == 0) {
    			result.Append(a);
    		}
    		for (let d = step;d < l; d += step) {
    			result.Append(a.Add(v.MulScalar(d/l)));
    		}
    		result.Append(b);
    	}
    	return result;
    }
    Filter(f) {
    	const result = new Paths();
    	let path = new Path();
    	for (const _ in this.v) {
    	    const v = this.v[_];
    		const fr = f.Filter(v);
    		if (fr.ok || (DRAW_HIDDEN_LINES && _%8 < 4)) { // show hidden lines
    			path.Append(fr.w);
    		} else {
    			if (path && path.v.length > 1) {
    				result.Append(path);
    			}
    			path = new Path();
    		}
    	}
    	if (path && path.v.length > 1) {
    		result.Append(path);
    	}
    	return result;
    }
    Simplify(threshold) {
    	if (this.v.length < 3) {
    		return this;
    	}
    	const a = this.v[0];
    	const b = this.v[this.v.length-1];
    	let index = -1;
    	let distance = 0.0;
    	for (let i = 1; i < this.v.length-1; i++) {
    		const d = this.v[i].SegmentDistance(a, b);
    		if (d > distance) {
    			index = i;
    			distance = d;
    		}
    	}
    	if (distance > threshold) {
    		const r1 = new Path(this.v.slice(0, index+1)).Simplify(threshold);
    		const r2 = new Path(this.v.slice(index)).Simplify(threshold);
    		return new Path([...r1.v.slice(0, r1.v.length-1), ...r2.v]);
    	} else {
    		return new Path([a, b]);
    	}
    }
}

class Paths {
    constructor(v = []) {
        this.v = v;
    }
    Append(v) {
        if (Array.isArray(v)) {
            this.v = this.v.concat(v);   
        } else {
            this.v.push(v);
        }
    }
    BoundingBox() {
    	let box = this.v[0].BoundingBox();
    	for (const _ in this.v) {
    	    const path = this.v[_];
    		box = box.Extend(path.BoundingBox());
    	}
    	return box;
    }
    Transform(matrix) {
    	const result = new Paths();
    	for (const _ in this.v) {
    	    const path = this.v[_];
    		result.Append(path.Transform(matrix));
    	}
    	return result;
    }
    Chop(step) {
    	const result = new Paths();
    	for (const _ in this.v) {
    	    const path = this.v[_];
    		result.Append(path.Chop(step));
    	}
    	return result;
    }
    Filter(f) {
    	const result = new Paths();
    	for (const _ in this.v) {
    	    const path = this.v[_];
    		result.Append(path.Filter(f).v);
    	}
    	return result;
    }
    Simplify(threshold) {
    	const result = new Paths();
    	for (const _ in this.v) {
    	    const path = this.v[_];
    		result.Append(path.Simplify(threshold));
    	}
    	return result;
    }
}
//
// https://github.com/fogleman/ln/blob/master/ln/box.go
//
class Box {
    constructor(min, max) {
    	this.Min = min;
    	this.Max = max;
    }
    Anchor(anchor) { return this.Min.Add(this.Size().Mul(anchor)); }
    Center() { return this.Anchor(new Vector(0.5, 0.5, 0.5)); }
    Size() { return this.Max.Sub(a.Min); }
    Contains(b) {
    	return this.Min.X <= b.X && this.Max.X >= b.X &&
    		this.Min.Y <= b.Y && this.Max.Y >= b.Y &&
    		this.Min.Z <= b.Z && this.Max.Z >= b.Z;
    }
    Extend(b) { return new Box(this.Min.Min(b.Min), this.Max.Max(b.Max)); }
    Intersect(r) {  	
        let x1 = (this.Min.X - r.Origin.X) / r.Direction.X;
    	let y1 = (this.Min.Y - r.Origin.Y) / r.Direction.Y;
    	let z1 = (this.Min.Z - r.Origin.Z) / r.Direction.Z;
    	let x2 = (this.Max.X - r.Origin.X) / r.Direction.X;
    	let y2 = (this.Max.Y - r.Origin.Y) / r.Direction.Y;
    	let z2 = (this.Max.Z - r.Origin.Z) / r.Direction.Z;
    	if (x1 > x2) {
    		const tmp = x1; x1 = x2, x2 = tmp;
    	}
    	if (y1 > y2) {
    		const tmp = y1; y1 = y2, y2 = tmp;
    	}
    	if (z1 > z2) {
    		const tmp = z1; z1 = z2, z2 = tmp;
    	}
    	const t1 = Math.max(Math.max(x1, y1), z1);
    	const t2 = Math.min(Math.min(x2, y2), z2);
    	return {t1, t2};
    }
    Partition(axis, point) {
        let left, right;
    	switch (axis) {
        	case AxisX:
        		left = this.Min.X <= point;
        		right = this.Max.X >= point;
        	case AxisY:
        		left = this.Min.Y <= point;
        		right = this.Max.Y >= point;
        	case AxisZ:
        		left = this.Min.Z <= point;
        		right = this.Max.Z >= point;
    	}
    	return {left, right};
    }
}
function BoxForShapes(shapes) {
	if (shapes.length == 0) {
		return new Box();
	}
	let box = shapes[0].BoundingBox();
	for (const _ in shapes) {
	    const shape = shapes[_];
		box = box.Extend(shape.BoundingBox());
	}
	return box;
}
function BoxForTriangles(shapes) {
	if (shapes.length == 0) {
		return new Box();
	}
	let box = shapes[0].BoundingBox();
	for (const _ in shapes) {
	    const shape = shapes[_];
		box = box.Extend(shape.BoundingBox());
	}
	return box;
}
function BoxForVectors(vectors) {
	if (vectors.length == 0) {
		return new Box();
	}
	let min = vectors[0]
	let max = vectors[0]
	for (const _ in vectors) {
	    const v = vectors[_];
		min = min.Min(v);
		max = max.Max(v);
	}
	return new Box(min, max);
}
//
// https://github.com/fogleman/ln/blob/master/ln/box.go
//
class Hit {
    constructor(shape, T) {
        this.Shape = shape;
        this.T = T;
    }
    Ok() {
    	return this.T < INF;
    }
    Min(b) {
    	if (this.T <= b.T) {
    		return this;
    	}
    	return b;
    }
    Max(b) {
    	if (this.T > b.T) {
    		return a;
    	}
    	return b;
    }
}
const NoHit = new Hit(false, INF);
//
// https://github.com/fogleman/ln/blob/master/ln/filter.go
//
const ClipBox = new Box(new Vector(-1, -1, -1), new Vector(1, 1, 1));
class ClipFilter {
    constructor(m, eye, scene) {
        this.Matrix = m;
        this.Eye = eye;
        this.Scene = scene;
    }
    Filter(v) {
    	const w = this.Matrix.MulPositionW(v);
    	if (!this.Scene.Visible(this.Eye, v)) {
    		return {w, ok: false};
    	}
    	if (!ClipBox.Contains(w)) {
    		return {w, ok: false};
    	}
    	return {w, ok: true};
    }
}
//  
// https://github.com/fogleman/ln/blob/master/ln/triangle.go
//
class Triangle {
    constructor(v1, v2, v3) {
    	this.V1 = v1; this.V2 = v2; this.V3 = v3;
    	this.UpdateBoundingBox();
    }
    UpdateBoundingBox() {
    	const min = this.V1.Min(this.V2).Min(this.V3);
    	const max = this.V1.Max(this.V2).Max(this.V3);
    	this.Box = new Box(min, max);
    }
    Compile() {}
    BoundingBox() { return this.Box; }
    Contains(v, f) { return false; }
    Intersect(r) {
    	const e1 = this.V2.Sub(this.V1);
    	const e2 = this.V3.Sub(this.V1);
    	const p =r.Direction.Cross(e2);
    	const det = e1.Dot(p);
    	if (det > -EPS && det < EPS) {
    		return NoHit;
    	}
    	const inv = 1 / det;
    	const t = r.Origin.Sub(this.V1);
    	const u = t.Dot(p) * inv;
    	if (u < 0 || u > 1) {
    		return NoHit;
    	}
    	const q = t.Cross(e1);
    	const v = r.Direction.Dot(q) * inv;
    	if (v < 0 || u+v > 1) {
    		return NoHit;
    	}
    	const d = e2.Dot(q) * inv;
    	if (d < EPS) {
    		return NoHit;
    	}
    	return new Hit(t, d);
    }
    Paths() { return new Paths([
        new Path([this.V1, this.V2]), 
        new Path([this.V2, this.V3]), 
        new Path([this.V3, this.V1])]);
    }
}
//  
// https://github.com/fogleman/ln/blob/master/ln/cube.go
//
class Cube {
    constructor(min, max) {
        this.Min = min;
        this.Max = max;
        this.Box = new Box(min, max);
    }
    Compile() {}
    BoundingBox() {
    	return this.Box;
    }
    Contains(v, f) {
    	if (v.X < this.Min.X-f || v.X > this.Max.X+f) {
    		return false;
    	}
    	if (v.Y < this.Min.Y-f || v.Y > this.Max.Y+f) {
    		return false;
    	}
    	if (v.Z < this.Min.Z-f || v.Z > this.Max.Z+f) {
    		return false;
    	}
    	return true;
    }
    Intersect(r) {
    	let n = this.Min.Sub(r.Origin).Div(r.Direction);
    	let f = this.Max.Sub(r.Origin).Div(r.Direction);
    	const v = n.Min(f); f = n.Max(f);
    	const t0 = Math.max(Math.max(v.X, v.Y), v.Z);
    	const t1 = Math.min(Math.min(f.X, f.Y), f.Z);
    	if (t0 < 1e-3 && t1 > 1e-3) {
    		return new Hit(this, t1);
    	}
    	if (t0 >= 1e-3 && t0 < t1) {
    		return new Hit(this, t0);
    	}
    	return NoHit;
    }
    Paths() {
    	const x1 = this.Min.X, y1 = this.Min.Y, z1 = this.Min.Z;
    	const x2 = this.Max.X, y2 = this.Max.Y, z2 = this.Max.Z;
    	const paths = new Paths([
    		new Path([new Vector(x1, y1, z1), new Vector(x1, y1, z2)]),
     		new Path([new Vector(x1, y1, z1), new Vector(x1, y2, z1)]),
     		new Path([new Vector(x1, y1, z1), new Vector(x2, y1, z1)]),
     		new Path([new Vector(x1, y1, z2), new Vector(x1, y2, z2)]),
     		new Path([new Vector(x1, y1, z2), new Vector(x2, y1, z2)]),
     		new Path([new Vector(x1, y2, z1), new Vector(x1, y2, z2)]),
     		new Path([new Vector(x1, y2, z1), new Vector(x2, y2, z1)]),
     		new Path([new Vector(x1, y2, z2), new Vector(x2, y2, z2)]),
     		new Path([new Vector(x2, y1, z1), new Vector(x2, y1, z2)]),
     		new Path([new Vector(x2, y1, z1), new Vector(x2, y2, z1)]),
     		new Path([new Vector(x2, y1, z2), new Vector(x2, y2, z2)]),
     		new Path([new Vector(x2, y2, z1), new Vector(x2, y2, z2)]),
    	]);
    	return paths;
    }
}
//
// https://github.com/fogleman/ln/blob/master/examples/csg.go
//
const Intersection = 0;
const Difference = 1;

class BooleanShape {
    constructor(op, A, B) {
        this.Op = op;
        this.A = arguments[1];
        this.B = arguments[2];
        for (let i=3; i<arguments.length; i++) {
		    this.A = new BooleanShape(op, this.A, arguments[i]);
	    }
    }
    Compile() {}
    BoundingBox() {
    	// TODO: fix this
    	const a = this.A.BoundingBox();
    	const b = this.B.BoundingBox();
    	return a.Extend(b);
    }
    Contains(v, f) {
    	f = 1e-3;
    	switch (this.Op) {
        	case Intersection:
        		return this.A.Contains(v, f) && this.B.Contains(v, f);
        	case Difference:
        		return this.A.Contains(v, f) && !this.B.Contains(v, -f);
    	}
    	return false;
    }
    Intersect(r) {
    	const h1 = this.A.Intersect(r);
    	const h2 = this.B.Intersect(r);
    	const h = h1.Min(h2);
    	const v = r.Position(h.T);
    	if (!h.Ok() || this.Contains(v, 0)) {
    		return h;
    	}
    	return this.Intersect(new Ray(r.Position(h.T+0.01), r.Direction));
    }
    Paths() {
    	let p = this.A.Paths();
    	p.Append(this.B.Paths().v);
    	return p.Chop(0.01).Filter(this);
    }
    Filter(w) {
    	return {w, ok: this.Contains(w, 0)};
    }
}
//
// https://github.com/fogleman/ln/blob/master/ln/cylinder.go
//
class Cylinder {
    constructor(radius, z0, z1) {
	    this.Radius = radius;
	    this.Z0 = z0;
	    this.Z1 = z1;
    }
    Compile() {}
    BoundingBox() {
    	const r = this.Radius;
    	return new Box(new Vector(-r, -r, this.Z0), new Vector(r, r, this.Z1));
    }
    Contains(v, f) {
    	const xy = new Vector(v.X, v.Y, 0);
    	if (xy.Length() > this.Radius+f) {
    		return false;
    	}
    	return v.Z >= this.Z0-f && v.Z <= this.Z1+f;
    }
    Intersect(ray)  {
    	const r = this.Radius;
    	const o = ray.Origin;
    	const d = ray.Direction;
    	const a = d.X*d.X + d.Y*d.Y;
    	const b = 2*o.X*d.X + 2*o.Y*d.Y;
    	const c = o.X*o.X + o.Y*o.Y - r*r;
    	const q = b*b - 4*a*c;
    	if (q < 0) {
    		return NoHit;
    	}
    	const s = Math.sqrt(q);
    	let t0 = (-b + s) / (2 * a);
    	let t1 = (-b - s) / (2 * a);
    	if (t0 > t1) {
    	    const tmp = t0;
    	    t0 = t1; t1 = tmp;
    	}
    	const z0 = o.Z + t0*d.Z;
    	const z1 = o.Z + t1*d.Z;
    	if (t0 > 1e-6 && this.Z0 < z0 && z0 < this.Z1) {
    		return new Hit(this, t0);
    	}
    	if (t1 > 1e-6 && this.Z0 < z1 && z1 < this.Z1) {
    		return new Hit(this, t1);
    	}
    	return NoHit;
    
    }
    Paths()  {
    	const result = new Paths();
    	for (let a = 0; a < 360; a += 10) {
    		const x = this.Radius * Math.cos(Radians(a));
    		const y = this.Radius * Math.sin(Radians(a));
    		result.Append(new Path([new Vector(x, y, this.Z0), new Vector(x, y, this.Z1)]));
    	}
    	return result;
    }
}

class OutlineCylinder extends Cylinder {
    constructor(eye, up, radius, z0, z1) {
        super(radius, z0, z1);
        this.Eye = eye;
        this.Up = up;
    }
    Paths() {
    	let center = new Vector(0, 0, this.Z0);
    	let hyp = center.Sub(this.Eye).Length();
    	let opp = this.Radius;
    	let theta = Math.asin(opp / hyp);
    	let adj = opp / Math.tan(theta);
    	let d = Math.cos(theta) * adj;
    	// r := math.Sin(theta) * adj
    	let w = center.Sub(this.Eye).Normalize();
    	let u = w.Cross(this.Up).Normalize();
    	const c0 = this.Eye.Add(w.MulScalar(d));
    	const a0 = c0.Add(u.MulScalar(this.Radius * 1.01));
    	const b0 = c0.Add(u.MulScalar(-this.Radius * 1.01));
    
    	center = new Vector(0, 0, this.Z1);
    	hyp = center.Sub(this.Eye).Length();
    	opp = this.Radius;
    	theta = Math.asin(opp / hyp);
    	adj = opp / Math.tan(theta);
    	d = Math.cos(theta) * adj;
    	// r = math.Sin(theta) * adj
    	w = center.Sub(this.Eye).Normalize();
    	u = w.Cross(this.Up).Normalize();
    	const c1 = this.Eye.Add(w.MulScalar(d));
    	const a1 = c1.Add(u.MulScalar(this.Radius * 1.01));
    	const b1 = c1.Add(u.MulScalar(-this.Radius * 1.01));
    
    	const p0 = new Path(), p1 = new Path();
    	for (let a = 0; a < 360; a++) {
    		const x = this.Radius * Math.cos(Radians(a));
    		const y = this.Radius * Math.sin(Radians(a));
    		p0.Append(new Vector(x, y, this.Z0));
    		p1.Append(new Vector(x, y, this.Z1));
    	}
    	return new Paths([
    		p0,
    		p1,
    		new Path([new Vector(a0.X, a0.Y, this.Z0), new Vector(a1.X, a1.Y, this.Z1)]),
    		new Path([new Vector(b0.X, b0.Y, this.Z0), new Vector(b1.X, b1.Y, this.Z1)]),
    		]);
    }
}

function NewTransformedOutlineCylinder(eye, up, v0, v1, radius) {
	const d = v1.Sub(v0);
	const z = d.Length();
	const a = Math.acos(d.Normalize().Dot(up));
	let m = Translate(v0);
	if (a != 0) {
		const u = d.Cross(up).Normalize();
		m = Rotate(u, a).Translate(v0);
	}
	const c = new OutlineCylinder(m.Inverse().MulPosition(eye), up, radius, 0, z);
	return new TransformedShape(c, m);
}
//
// https://github.com/fogleman/ln/blob/master/ln/sphere.go
//
class Sphere {
    constructor(center, radius) {
        this.Center = center;
        this.Radius = radius;
        
        const min = new Vector(center.X - radius, center.Y - radius, center.Z - radius);
    	const max = new Vector(center.X + radius, center.Y + radius, center.Z + radius);
    	this.Box = new Box(min, max);
    }
    Compile() {}
    BoundingBox() {
        return this.Box;
    }
    Contains(v, f) {
    	return v.Sub(this.Center).Length() <= this.Radius+f;
    }
    Intersect(r) {
    	const radius = this.Radius;
    	const to = r.Origin.Sub(this.Center);
    	const b = to.Dot(r.Direction);
    	const c = to.Dot(to) - radius*radius;
    	let d = b*b - c;
    	if (d > 0) {
    		d = Math.sqrt(d);
    		const t1 = -b - d;
    		if (t1 > 1e-2) {
    			return new Hit(this, t1);
    		}
    		const t2 = -b + d;
    		if (t2 > 1e-2) {
    			return new Hit(this, t2);
    		}
    	}
    	return NoHit
    }
    Paths3() {
    	const paths = new Paths();
    	for (let i = 0; i < 20000; i++) {
    		let v = RandomUnitVector();
    		v = v.MulScalar(this.Radius).Add(this.Center);
    		paths.Append(new Path([v, v]));
    	}
    	return paths;
    }
    Paths2() {
    	var equator = new Path();
    	for (let lng = 0; lng <= 360; lng++) {
    		const v = LatLngToXYZ(0, lng, this.Radius);
    		equator.Append(v);
    	}
    	var paths = new Paths();
    	for (let i = 0; i < 100; i++) {
    		const m = Identity();
    		for (let j = 0; j < 3; j++) {
    			const v = RandomUnitVector();
    			m = m.Rotate(v, Math.random()*2*Math.PI);
    		}
    		m = m.Translate(this.Center);
    		paths.Append(equator.Transform(m));
    	}
    	return paths;
    }
    Paths() {
    	const paths = new Paths();
    	const n = 10;
    	const o = 10;
    	for (let lat = -90 + o; lat <= 90-o; lat += n) {
    		const path = new Path();
    		for (let lng = 0; lng <= 360; lng++) {
    			const v = LatLngToXYZ(lat, lng, this.Radius).Add(this.Center);
    			path.Append(v);
    		}
    		paths.Append(path);
    	}
    	for (let lng = 0; lng < 360; lng += n) {
    		const path = new Path();
    		for (let lat = -90 + o; lat <= 90-o; lat++) {
    			const v = LatLngToXYZ(lat, lng, this.Radius).Add(this.Center);
    			path.Append(v);
    		}
    		paths.Append(path);
    	}
    	return paths;
    }

}

function LatLngToXYZ(lat, lng, radius) {
	const latr = Radians(lat);
	const lngr = Radians(lng);
	return new Vector(  radius * Math.cos(latr) * Math.cos(lngr),
	                    radius * Math.cos(latr) * Math.sin(lngr),
	                    radius * Math.sin(latr) );
}

class OutlineSphere extends Sphere {
    constructor(eye, up, center, radius) {
        super(center, radius);
        this.Eye = eye;
	    this.Up = up;
    }
    Paths() {
    	const path = new Path();
    	const center = this.Center;
    	const radius = this.Radius;
    
    	const hyp = center.Sub(this.Eye).Length();
    	const opp = radius;
    	const theta = Math.asin(opp / hyp);
    	const adj = opp / Math.tan(theta);
    	const d = Math.cos(theta) * adj;
    	const r = Math.sin(theta) * adj;
    
    	const w = center.Sub(this.Eye).Normalize()
    	const u = w.Cross(this.Up).Normalize()
    	const v = w.Cross(u).Normalize()
    	const c = this.Eye.Add(w.MulScalar(d))
    	for (let i = 0; i <= 360; i++) {
    		const a = Radians(i);
    		let p = c;
    		p = p.Add(u.MulScalar(Math.cos(a) * r));
    		p = p.Add(v.MulScalar(Math.sin(a) * r));
    		path.Append(p);
    	}
    	return new Paths([path]);
    }
}
//
// https://github.com/fogleman/ln/blob/master/ln/function.go
//
const Above = 0;
const Below = 1;

class FunctionShape {
    constructor(func, box, direction) {
	    this.Func = func;
	    this.Box = box;
	    this.Direction = direction;
    }
    Compile() {}
    BoundingBox() { return this.Box; }
    Contains(v, eps) {
    	if (this.Direction == Below) {
    		return v.Z < this.Func(v.X, v.Y);
    	} else {
    		return v.Z > this.Func(v.X, v.Y);
    	}
    }
    Intersect(ray) {
    	const step = 1.0 / 64;
    	const sign = this.Contains(ray.Position(step), 0);
    	for (let t = step; t < 10; t += step) {
    		const v = ray.Position(t);
    		if (this.Contains(v, 0) != sign && this.Box.Contains(v)) {
    			return new Hit(this, t);
    		}
    	}
    	return NoHit
    }
    Paths3() {
    	const path = new Path();
    	const n = 10000
    	for (let i = 0; i < n; i++) {
    		const t = i / n;
    		const r = 8 - Math.pow(t, 0.1)*8;
    		const x = Math.cos(Radians(t*2*Math.PI*3000)) * r;
    		const y = Math.sin(Radians(t*2*Math.PI*3000)) * r;
    		let z = this.Func(x, y);
    		z = Math.min(z, this.Box.Max.Z);
    		z = Math.max(z, this.Box.Min.Z);
    		path.Append(new Vector(x, y, z));
    	}
    	// return append(f.Paths2(), path)
    	return new Paths([path]);
    }
    Paths() {
    	const paths = new Paths();
    	const fine = 1.0 / 256;
    	for (let a = 0; a < 360; a += 5) {
    		const path = new Path();
    		for (let r = 0.0; r <= 8.0; r += fine) {
    			let x = Math.cos(Radians(a)) * r;
    			let y = Math.sin(Radians(a)) * r;
    			let z = this.Func(x, y);
    			const o = -Math.pow(-z, 1.4);
    			x = Math.cos(Radians(a)-o) * r;
    			y = Math.sin(Radians(a)-o) * r;
    			z = Math.min(z, this.Box.Max.Z);
    			z = Math.max(z, this.Box.Min.Z);
    			path.Append(new Vector(x, y, z));
    		}
    		paths.Append(path);
    	}
    	return paths;
    }
    Paths1() {
    	const paths = new Paths();
    	const step = 1.0 / 8;
    	const fine = 1.0 / 64;
    	for (let x = this.Box.Min.X; x <= this.Box.Max.X; x += step) {
    		const path = new Path();
    		for (let y = this.Box.Min.Y; y <= this.Box.Max.Y; y += fine) {
    			let z = this.Func(x, y);
    			z = Math.min(z, this.Box.Max.Z);
    			z = Math.max(z, this.Box.Min.Z);
    			path.Append(new Vector(x, y, z));
    		}
    		paths.Append(path);
    	}
    	for (let y = this.Box.Min.Y; y <= this.Box.Max.Y; y += step) {
    		const path = new Path();
    		for (let x = this.Box.Min.X; x <= this.Box.Max.X; x += fine) {
    			let z = this.Func(x, y);
    			z = Math.min(z, this.Box.Max.Z);
    			z = Math.max(z, this.Box.Min.Z);
    			path.Append(new Vector(x, y, z));
    		}
    		paths.Append([path]);
    	}
    	return paths
    }
}
