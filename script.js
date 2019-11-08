// https://github.com/fogleman/ln
// https://turtletoy.net/turtle/9d2e6913b4
let paths, scene;

const width = 300*2;
const height = 210*2;
const DRAW_HIDDEN_LINES = false;

let eye,center,up,fovy,znear=0.1,zfar=100,step=0.01;


// --------------------------------------------------------
function setup()
{

	createCanvas(width,height,SVG);
	createScene("cosinus-city");
}

// --------------------------------------------------------
function draw()
{
	renderScene();
  	noLoop();
	save();
}

// --------------------------------------------------------
function renderScene()
{
	paths = scene.Render(eye, center, up, width, height, fovy, znear, zfar, step);
    if (paths) 
    {
    	noFill();
    	stroke(0);
    	translate(width/2, height/2);

        for (let k=0; k<paths.v.length; k++) 
        {
            const path = paths.v[k];
            beginShape();
            for (let j=0; j<path.v.length; j++)
                vertex(path.v[j].X, path.v[j].Y);
            endShape();
        }
  	}
}

// --------------------------------------------------------
function createScene(name)
{
	scene = new Scene();

	if (name == "spheres")
	{
		eye = new Vector(8, 8, 8);
		center = new Vector(0, 0, 1);
		up = new Vector(0, 0, 1);
		fovy = 20.0;
		const n = 3;
    	for (let x = -n; x <= n; x++) {
    		for (let y = -n; y <= n; y++) {
    			const z = Math.random() * 3;
    			const v = new Vector(x, y, z);
    			// const sphere = new Cube(v, new Vector(v.X+.4,v.Y+.4,v.Z+.4));// .45);
    			const sphere = new Sphere(v, .45);
    			// const sphere = new OutlineSphere(eye, up, v, .45);
    			scene.Add(sphere);
    		}
    	}

	}
	else if (name == "skyscrapers")
	{

		eye = new Vector(1.75, 1.25, 6);
		center = new Vector(0, 0, 0);
		up = new Vector(0, 0, 1);
		fovy = 100.0;

		let p = Math.random()*0.25 + 0.2;				
		let shape = new Cube(new Vector(0 - p, 0 - p, 0), new Vector(0 + p, 0 + p, 2));

		let n = 7;
		for (let x=-n; x<=n ; x++)
		{
			for (let y=-n; y<=n ; y++)
			{
				let p = Math.random()*0.25 + 0.2;				
				let dx = Math.random()*0.5 - 0.25;				
				let dy = Math.random()*0.5 - 0.25;				
				let fx = (x) + dx*0;
    			let fy = (y) + dy*0;
				let fz = Math.random()*3 + 1;

    			let shape = new Cube(new Vector(fx - p, fy - p, 0), new Vector(fx + p, fy + p, fz));
    			if (!(x == 2 && y == 1)) {
    			    scene.Add(shape);
    			}
			}
		}
	}
	else if (name == "cosinus-city")
	{
		eye = new Vector(0.01, 0, 10);
		center = new Vector(0, 0, 0);
		up = new Vector(0, 0, 1);
		fovy = 100.0;
		let angle = 0;
		let r = 7.5;
		while (r>=0)
		{
		 	let x = r*cos( Radians(angle) );
   			let y = r*sin( Radians(angle) );
   			let sx = 0.25;
			let sy = 0.25;				
			let sz = Math.random()*2.5+1;
			
			let shape = new Cube(new Vector(x - sx, y - sy, 0), new Vector(x + sx, y + sy, sz));
			scene.Add(shape);
			r = r - 0.010;
   			angle += 2.5;
   		}

	}
}

