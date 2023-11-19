use std::f32;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub enum Fields {
    UField,
    VField,
    SField,
}

#[wasm_bindgen]
pub struct FluidSimulator {
    density: f32,
    gravity: f32,
    over_relaxation: f32,
    num_x: usize,
    num_y: usize,
    num_cells: usize,
    h: f32,
    u: Vec<f32>,
    m: Vec<f32>,
    s: Vec<f32>,
    v: Vec<f32>,
    new_u: Vec<f32>,
    new_v: Vec<f32>,
    new_m: Vec<f32>,
    p: Vec<f32>,
    drag: Vec<f32>,
    lift: Vec<f32>,
}

#[wasm_bindgen]
impl FluidSimulator {
    #[wasm_bindgen(constructor)]
    pub fn new(density: f32, num_x: usize, num_y: usize, h: f32) -> Self {
        let num_x = num_x + 2;
        let num_y = num_y + 2;
        let num_cells = num_x * num_y;
        let u = vec![0.0; num_cells];
        let v = vec![0.0; num_cells];
        let new_u = vec![0.0; num_cells];
        let new_v = vec![0.0; num_cells];
        let p = vec![0.0; num_cells];
        let s = vec![0.0; num_cells];
        let m = vec![1.0; num_cells];
        let new_m = vec![0.0; num_cells];
        let drag = vec![0.0; num_cells];
        let lift = vec![0.0; num_cells];
        FluidSimulator {
            density,
            gravity: -9.81,
            over_relaxation: 1.90,
            num_x,
            num_y,
            num_cells,
            h,
            u,
            m,
            s,
            v,
            new_u,
            new_v,
            new_m,
            p,
            drag,
            lift,
        }
    }

    pub fn integrate(&mut self, dt: f32) {
        let n = self.num_y;
        for i in 1..self.num_x {
            for j in 1..self.num_y - 1 {
                if self.s[i * n + j] != 0.0 && self.s[i * n + j - 1] != 0.0 {
                    self.v[i * n + j] += self.gravity * dt;
                }
            }
        }
    }

    pub fn solve_incompressibility(&mut self, num_iters: usize, dt: f32) {
        let n = self.num_y;
        let cp = self.density * self.h / dt;

        for _ in 0..num_iters {
            for i in 1..self.num_x - 1 {
                for j in 1..self.num_y - 1 {
                    if self.s[i * n + j] == 0.0 {
                        continue;
                    }

                    let sx0 = self.s[(i - 1) * n + j];
                    let sx1 = self.s[(i + 1) * n + j];
                    let sy0 = self.s[i * n + j - 1];
                    let sy1 = self.s[i * n + j + 1];
                    let s = sx0 + sx1 + sy0 + sy1;
                    if s == 0.0 {
                        continue;
                    }

                    let div = self.u[(i + 1) * n + j] - self.u[i * n + j] + self.v[i * n + j + 1]
                        - self.v[i * n + j];

                    let mut p = -div / s;
                    p *= self.over_relaxation;
                    self.p[i * n + j] += cp * p;

                    self.u[i * n + j] -= sx0 * p;
                    self.u[(i + 1) * n + j] += sx1 * p;
                    self.v[i * n + j] -= sy0 * p;
                    self.v[i * n + j + 1] += sy1 * p;
                }
            }
        }
    }

    pub fn extrapolate(&mut self) {
        let n = self.num_y;
        for i in 0..self.num_x {
            self.u[i * n] = self.u[i * n + 1];
            self.u[i * n + n - 1] = self.u[i * n + n - 2];
        }
        for j in 0..self.num_y {
            self.v[j] = self.v[n + j];
            self.v[(self.num_x - 1) * n + j] = self.v[(self.num_x - 2) * n + j];
        }
    }

    pub fn sample_field(&self, x: f32, y: f32, field: Fields) -> f32 {
        let n = self.num_y;
        let h = self.h;
        let h1 = 1.0 / h;
        let h2 = 0.5 * h;

        let x = x.max(h).min(self.num_x as f32 * h);
        let y = y.max(h).min(self.num_y as f32 * h);

        let mut dx = 0.0;
        let mut dy = 0.0;

        let f: &[f32];

        match field {
            Fields::UField => {
                f = &self.u;
                dy = h2;
            }
            Fields::VField => {
                f = &self.v;
                dx = h2;
            }
            Fields::SField => {
                f = &self.m;
                dx = h2;
                dy = h2;
            }
        }

        let x0 = (x - dx) * h1;
        let tx = (x - dx - x0 * h) * h1;
        let x0 = x0.min(self.num_x as f32 - 1.0) as usize;
        let x1 = x0 + 1;

        let y0 = (y - dy) * h1;
        let ty = (y - dy - y0 * h) * h1;
        let y0 = y0.min(self.num_y as f32 - 1.0) as usize;
        let y1 = y0 + 1;

        let sx = 1.0 - tx;
        let sy = 1.0 - ty;

        sx * sy * f[x0 * n + y0]
            + tx * sy * f[x1 * n + y0]
            + tx * ty * f[x1 * n + y1]
            + sx * ty * f[x0 * n + y1]
    }

    pub fn average_field(&self, i: usize, j: usize, field: Fields) -> f32 {
        let n = self.num_y;
        match field {
            Fields::UField => {
                (self.u[i * n + j - 1]
                    + self.u[i * n + j]
                    + self.u[(i + 1) * n + j - 1]
                    + self.u[(i + 1) * n + j])
                    * 0.25
            }
            Fields::VField => {
                (self.v[(i - 1) * n + j]
                    + self.v[i * n + j]
                    + self.v[(i - 1) * n + j + 1]
                    + self.v[i * n + j + 1])
                    * 0.25
            }
            Fields::SField => 0.0,
        }
    }

    pub fn advect_vel(&mut self, dt: f32) {
        self.new_u.copy_from_slice(&self.u);
        self.new_v.copy_from_slice(&self.v);

        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;

        for i in 1..self.num_x {
            for j in 1..self.num_y {
                // u component
                if self.s[i * n + j] != 0.0 && self.s[(i - 1) * n + j] != 0.0 && j < self.num_y - 1
                {
                    let x = i as f32 * h;
                    let y = j as f32 * h + h2;
                    let u = self.u[i * n + j];
                    let v = self.average_field(i, j, Fields::VField);
                    let x = x - dt * u;
                    let y = y - dt * v;
                    let u = self.sample_field(x, y, Fields::UField);
                    self.new_u[i * n + j] = u;
                }
                // v component
                if self.s[i * n + j] != 0.0 && self.s[i * n + j - 1] != 0.0 && i < self.num_x - 1 {
                    let x = i as f32 * h + h2;
                    let y = j as f32 * h;
                    let u = self.average_field(i, j, Fields::UField);
                    let v = self.v[i * n + j];
                    let x = x - dt * u;
                    let y = y - dt * v;
                    let v = self.sample_field(x, y, Fields::VField);
                    self.new_v[i * n + j] = v;
                }
            }
        }

        self.u.copy_from_slice(&self.new_u);
        self.v.copy_from_slice(&self.new_v);
    }

    pub fn advect_smoke(&mut self, dt: f32) {
        self.new_m.copy_from_slice(&self.m);

        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;

        for i in 1..self.num_x - 1 {
            for j in 1..self.num_y - 1 {
                if self.s[i * n + j] != 0.0 {
                    let u = (self.u[i * n + j] + self.u[(i + 1) * n + j]) * 0.5;
                    let v = (self.v[i * n + j] + self.v[i * n + j + 1]) * 0.5;
                    let x = i as f32 * h + h2 - dt * u;
                    let y = j as f32 * h + h2 - dt * v;

                    self.new_m[i * n + j] = self.sample_field(x, y, Fields::SField);
                }
            }
        }
        self.m.copy_from_slice(&self.new_m);
    }

    //Lift and drag forces
    pub fn calculate_forces(&mut self) {
        let n = self.num_y;
        for i in 1..self.num_x - 1 {
            for j in 1..self.num_y - 1 {
                if self.s[i * n + j] == 0.0 {
                    continue;
                }

                // Calculate the velocity of the fluid at the current cell
                let velocity_x = self.u[i * n + j];
                let velocity_y = self.v[i * n + j];
                let velocity_magnitude =
                    f32::sqrt(f32::powf(velocity_x, 2.0) + f32::powf(velocity_y, 2.0));

                // Calculate the unit normal vector to the surface of the airfoil at the current cell
                let surface_normal_x = -self.u[i * n + j] / velocity_magnitude;
                let surface_normal_y = -self.v[i * n + j] / velocity_magnitude;

                // Calculate the lift force acting on the fluid at the current cell
                self.lift[i * n + j] =
                    self.density * f32::powf(velocity_magnitude, 2.0) * self.h * surface_normal_y;

                // Calculate the drag force acting on the fluid at the current cell
                self.drag[i * n + j] =
                    self.density * f32::powf(velocity_magnitude, 2.0) * self.h * surface_normal_x;
            }
        }
    }

    pub fn simulate(&mut self, dt: f32, num_iters: usize) {
        self.integrate(dt);

        self.p.iter_mut().for_each(|p| *p = 0.0);
        self.solve_incompressibility(num_iters, dt);

        self.extrapolate();
        self.advect_vel(dt);
        self.advect_smoke(dt);
        self.calculate_forces();
    }

    pub fn set_gravity(&mut self, g: f32) {
        self.gravity = g;
    }

    pub fn get_h(&mut self) -> f32 {
        self.h
    }

    pub fn get_num_x(&mut self) -> usize {
        self.num_x
    }

    pub fn get_num_y(&mut self) -> usize {
        self.num_y
    }

    pub fn get_num_cells(&mut self) -> usize {
        self.num_cells
    }

    pub fn get_p(&mut self, i: usize) -> f32 {
        self.p[i]
    }

    pub fn get_v(&mut self, i: usize) -> f32 {
        self.v[i]
    }

    pub fn set_p(&mut self, index:usize, value:f32) {
        self.p[index] = value;
    }

    pub fn set_h(&mut self, h: f32) {
        self.h = h;
    }

    pub fn set_u(&mut self, index:usize, value:f32) {
        self.u[index] = value;
    }

    pub fn set_m(&mut self, index:usize, value:f32) {
        self.m[index] = value;
    }

    pub fn set_s(&mut self, index:usize, value:f32) {
        self.s[index] = value;
    }

    pub fn set_v(&mut self, index:usize, value:f32) {
        self.v[index] = value;
    }

    pub fn get_m(&mut self, index:usize)->f32{
        self.m[index]
    }

    pub fn get_s(&mut self, index:usize)->f32{
        self.s[index]
    }

    pub fn get_u(&mut self, index:usize)->f32{
        self.u[index]
    }


    //add all getter and setters
}
