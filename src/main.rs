#![allow(dead_code)]
#![allow(non_snake_case)]
use geographiclib_rs::{DirectGeodesic, Geodesic, capability as caps};
#[allow(unused)]
// ============================================================
// curvas - Traducción (aproximada) del programa Python + curva geodésica + loxodrómica
use libm::{hypot, remainder};
use shapefile::{
    Point, Polyline, Writer,
    dbase::{FieldName, FieldValue, Record, TableWriterBuilder},
};
use std::f64::consts::PI;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::process;

// ============================================================
// CONSTANTES GLOBALES
//const DEBUG: bool = false;
const MAX_ITER: u32 = 200_000;

// Elipsoide GRS80
const GRS80_A: f64 = 6378137.0;
const GRS80_F: f64 = 1.0 / 298.2572221008827;

const EPSILON: f64 = 1E-15;
const PI_2: f64 = PI / 2.0;
const PI_4: f64 = PI / 4.0;

const MIN_FLOAT: f64 = f64::MIN_POSITIVE;
const MAX_FLOAT: f64 = f64::MAX;

const SQRT_MAX_FLOAT: f64 = 1.3407807929942596e154; // sqrt(MAX_FLOAT)
const SQRT_MIN_FLOAT: f64 = 1.4916681462400413e-154; // sqrt(MIN_FLOAT)

const EPSILON_MAQ: f64 = 1E-16;
const EPSILON_ANG: f64 = 1E-10;
const MAXSTEP_ANG: f64 = 0.08726646259971647; // 5 grados en radianes
const MAXSTEPRK45_ANG_INITIALVALUE: f64 = 0.0017453292519943296; // 0.1 grados en radianes
const MAXSTEP_H: f64 = 1E-5;
const MAXDIST: f64 = 10.0;

// ============================================================
// ESTRUCTURAS DE DATOS
#[derive(Debug, Clone)]
struct Pelipsoide {
    a: f64,
    f: f64,
    b: f64,
    c: f64,
    n: f64,
    n2: f64,
    n3: f64,
    n4: f64,
    n5: f64,
    n6: f64,
    e2: f64,
    e: f64,
    e4: f64,
    e6: f64,
    e8: f64,
    e10: f64,
    ep2: f64,
    ep: f64,
    r1: f64,
    r2: f64,
    r3: f64,
    rm: f64,
    rmu: f64,
}

#[derive(Debug, Clone)]
struct Platn6 {
    geod2rect: [f64; 7],
    rect2geod: [f64; 7],
    isom2geod: [f64; 7],
    geod2conf: [f64; 7],
    conf2geod: [f64; 7],
    geod2auth: [f64; 7],
    auth2geod: [f64; 7],
    geod2geoc: [f64; 7],
    geoc2geod: [f64; 7],
    geod2para: [f64; 7],
    para2geod: [f64; 7],
    conf2rect: [f64; 7],
    auth2rect: [f64; 7],
    rect2conf: [f64; 7],
    rect2auth: [f64; 7],
    geoc2rect: [f64; 7],
    rect2geoc: [f64; 7],
    para2rect: [f64; 7],
    rect2para: [f64; 7],
    auth2conf: [f64; 7],
    conf2auth: [f64; 7],
    geoc2conf: [f64; 7],
    conf2geoc: [f64; 7],
    para2conf: [f64; 7],
    conf2para: [f64; 7],
    geoc2auth: [f64; 7],
    auth2geoc: [f64; 7],
    para2auth: [f64; 7],
    auth2para: [f64; 7],
    para2geoc: [f64; 7],
    geoc2para: [f64; 7],
}

#[derive(Debug, Clone)]
struct PRK45 {
    tol: f64,
    h_max: f64,
    h_min: f64,
}

// ============================================================
// FUNCIONES MATEMÁTICAS MISC, MÁS ADELANTE APARECE LAM
fn tosemicirc(num: f64) -> f64 {
    let result = remainder(num, 2.0 * PI);
    if result == -PI { PI } else { result }
}

fn toquadcirc(num: f64) -> f64 {
    let mut num = tosemicirc(num);
    let sgn_num = sgn(num);
    num = num.abs();
    let result;
    if num > PI_2 {
        result = sgn_num * (PI - num);
    } else {
        result = sgn_num * num;
    }
    result
}

fn sincos(x: f64) -> (f64, f64) {
    x.sin_cos()
}

fn arccos(x: f64) -> f64 {
    let x_clipped = x.clamp(-1.0, 1.0);
    arg(x_clipped, (1.0 - x_clipped * x_clipped).sqrt())
}

fn arcsen(x: f64) -> f64 {
    let x_clipped = x.clamp(-1.0, 1.0);
    arg((1.0 - x_clipped * x_clipped).sqrt(), x_clipped)
}

fn sgn(num: f64) -> f64 {
    if num < 0.0 { -1.0 } else { 1.0 }
}

fn div(y: f64, x: f64) -> f64 {
    let sx = sgn(x);
    y / (x + sx * MIN_FLOAT * 2.0)
}

fn rad2deg(alpha: f64) -> f64 {
    alpha * 180.0 / PI
}

fn deg2rad(alpha: f64) -> f64 {
    alpha / 180.0 * PI
}

fn sqr(x: f64) -> f64 {
    if x.abs() < SQRT_MAX_FLOAT { if x.abs() > SQRT_MIN_FLOAT { x * x } else { MIN_FLOAT } } else { MAX_FLOAT }
}

fn arg(x: f64, y: f64) -> f64 {
    if y.abs() > EPSILON {
        let rho = x.hypot(y);
        2.0 * (div(y, rho + x)).atan()
    } else if x >= 0.0 {
        0.0
    } else {
        PI
    }
}

fn lam(chi: f64) -> f64 {
    let sgnchi = sgn(chi);
    let chi: f64 = chi.abs();
    let q: f64;
    if chi < PI_2 {
        q = (PI_4 + 0.5 * chi).tan().ln() * sgnchi;
    } else {
        q = 40.0 * sgnchi;
    }
    q
}

// ============================================================
// INICIALIZACIÓN DEL ELIPSOIDE
impl Pelipsoide {
    fn new(in_f: f64, in_a: f64) -> Self {
        let b = in_a * (1.0 - in_f);
        let c = in_a / (1.0 - in_f);
        let n = in_f / (2.0 - in_f);
        let n2 = n * n;
        let n3 = n2 * n;
        let n4 = n2 * n2;
        let n5 = n4 * n;
        let n6 = n5 * n;
        let e2 = in_f * (2.0 - in_f);
        let e = e2.sqrt();
        let e4 = e2 * e2;
        let e6 = e4 * e2;
        let e8 = e4 * e4;
        let e10 = e6 * e4;
        let ep2 = e2 / (1.0 - e2);
        let ep = ep2.sqrt();
        let r1 = (2.0 * in_a + b) / 3.0;
        let r2 = in_a * (0.5 + 0.5 * (div(1.0, e) - e) * e.atanh()).sqrt();
        let r3 = (in_a * in_a * b).cbrt();
        let rm = in_a / (1.0 + n);
        let rmu = rm * (1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0);
        Pelipsoide {
            a: in_a,
            f: in_f,
            b,
            c,
            n,
            n2,
            n3,
            n4,
            n5,
            n6,
            e2,
            e,
            e4,
            e6,
            e8,
            e10,
            ep2,
            ep,
            r1,
            r2,
            r3,
            rm,
            rmu,
        }
    }
}

// ============================================================
// INICIALIZACIÓN DE PLATN6
impl Platn6 {
    fn new(n: f64) -> Self {
        let n2 = n * n;
        let n3 = n2 * n;
        let n4 = n2 * n2;
        let n5 = n3 * n2;
        let n6 = n3 * n3;

        let mut p = Platn6 {
            geod2rect: [0.0; 7],
            rect2geod: [0.0; 7],
            isom2geod: [0.0; 7],
            geod2conf: [0.0; 7],
            conf2geod: [0.0; 7],
            geod2auth: [0.0; 7],
            auth2geod: [0.0; 7],
            geod2geoc: [0.0; 7],
            geoc2geod: [0.0; 7],
            geod2para: [0.0; 7],
            para2geod: [0.0; 7],
            conf2rect: [0.0; 7],
            auth2rect: [0.0; 7],
            rect2conf: [0.0; 7],
            rect2auth: [0.0; 7],
            geoc2rect: [0.0; 7],
            rect2geoc: [0.0; 7],
            para2rect: [0.0; 7],
            rect2para: [0.0; 7],
            auth2conf: [0.0; 7],
            conf2auth: [0.0; 7],
            geoc2conf: [0.0; 7],
            conf2geoc: [0.0; 7],
            para2conf: [0.0; 7],
            conf2para: [0.0; 7],
            geoc2auth: [0.0; 7],
            auth2geoc: [0.0; 7],
            para2auth: [0.0; 7],
            auth2para: [0.0; 7],
            para2geoc: [0.0; 7],
            geoc2para: [0.0; 7],
        };

        // GEOD2RECT
        p.geod2rect[1] = -3.0 / 2.0 * n + 9.0 / 16.0 * n3 - 3.0 / 32.0 * n5;
        p.geod2rect[2] = 15.0 / 16.0 * n2 - 15.0 / 32.0 * n4 + 135.0 / 2048.0 * n6;
        p.geod2rect[3] = -35.0 / 48.0 * n3 + 105.0 / 256.0 * n5;
        p.geod2rect[4] = 315.0 / 512.0 * n4 - 189.0 / 512.0 * n6;
        p.geod2rect[5] = -693.0 / 1280.0 * n5;
        p.geod2rect[6] = 1001.0 / 2048.0 * n6;

        // RECT2GEOD
        p.rect2geod[1] = 3.0 / 2.0 * n - 27.0 / 32.0 * n3 + 269.0 / 512.0 * n5;
        p.rect2geod[2] = 21.0 / 16.0 * n2 - 55.0 / 32.0 * n4 + 6759.0 / 4096.0 * n6;
        p.rect2geod[3] = 151.0 / 96.0 * n3 - 417.0 / 128.0 * n5;
        p.rect2geod[4] = 1097.0 / 512.0 * n4 - 15543.0 / 2560.0 * n6;
        p.rect2geod[5] = 8011.0 / 2560.0 * n5;
        p.rect2geod[6] = 293393.0 / 61440.0 * n6;

        // ISOM2GEOD
        p.isom2geod[1] = 4.0 * n - 32.0 / 3.0 * n2 + 124.0 / 5.0 * n3 - 3296.0 / 63.0 * n4 + 32476.0 / 315.0 * n5 - 30081056.0 / 155925.0 * n6;
        p.isom2geod[2] = 56.0 / 3.0 * n2 - 1984.0 / 15.0 * n3 + 65872.0 / 105.0 * n4 - 764096.0 / 315.0 * n5 + 85344488.0 / 10395.0 * n6;
        p.isom2geod[3] = 1792.0 / 15.0 * n3 - 149984.0 / 105.0 * n4 + 1085824.0 / 105.0 * n5 - 1801378688.0 / 31185.0 * n6;
        p.isom2geod[4] = 273856.0 / 315.0 * n4 - 931328.0 / 63.0 * n5 + 4510583296.0 / 31185.0 * n6;
        p.isom2geod[5] = 2137088.0 / 315.0 * n5 - 57822208.0 / 385.0 * n6;
        p.isom2geod[6] = 1232232448.0 / 22275.0 * n6;

        // GEOD2CONF
        p.geod2conf[1] = -2.0 * n + 2.0 / 3.0 * n2 + 4.0 / 3.0 * n3 - 82.0 / 45.0 * n4 + 32.0 / 45.0 * n5 + 4642.0 / 4725.0 * n6;
        p.geod2conf[2] = 5.0 / 3.0 * n2 - 16.0 / 15.0 * n3 - 13.0 / 9.0 * n4 + 904.0 / 315.0 * n5 - 1522.0 / 945.0 * n6;
        p.geod2conf[3] = -26.0 / 15.0 * n3 + 34.0 / 21.0 * n4 + 8.0 / 5.0 * n5 - 12686.0 / 2835.0 * n6;
        p.geod2conf[4] = 1237.0 / 630.0 * n4 - 12.0 / 5.0 * n5 - 24832.0 / 14175.0 * n6;
        p.geod2conf[5] = -734.0 / 315.0 * n5 + 109598.0 / 31185.0 * n6;
        p.geod2conf[6] = 444337.0 / 155925.0 * n6;

        // CONF2GEOD
        p.conf2geod[1] = 2.0 * n - 2.0 / 3.0 * n2 - 2.0 * n3 + 116.0 / 45.0 * n4 + 26.0 / 45.0 * n5 - 2854.0 / 675.0 * n6;
        p.conf2geod[2] = 7.0 / 3.0 * n2 - 8.0 / 5.0 * n3 - 227.0 / 45.0 * n4 + 2704.0 / 315.0 * n5 + 2323.0 / 945.0 * n6;
        p.conf2geod[3] = 56.0 / 15.0 * n3 - 136.0 / 35.0 * n4 - 1262.0 / 105.0 * n5 + 73814.0 / 2835.0 * n6;
        p.conf2geod[4] = 4279.0 / 630.0 * n4 - 332.0 / 35.0 * n5 - 399572.0 / 14175.0 * n6;
        p.conf2geod[5] = 4174.0 / 315.0 * n5 - 144838.0 / 6237.0 * n6;
        p.conf2geod[6] = 601676.0 / 22275.0 * n6;

        // GEOD2AUTH
        p.geod2auth[1] = -4.0 / 3.0 * n - 4.0 / 45.0 * n2 + 88.0 / 315.0 * n3 + 538.0 / 4725.0 * n4 + 20824.0 / 467775.0 * n5 - 44732.0 / 2837835.0 * n6;
        p.geod2auth[2] = 34.0 / 45.0 * n2 + 8.0 / 105.0 * n3 - 2482.0 / 14175.0 * n4 - 37192.0 / 467775.0 * n5 - 12467764.0 / 212837625.0 * n6;
        p.geod2auth[3] = -1532.0 / 2835.0 * n3 - 898.0 / 14175.0 * n4 + 54968.0 / 467775.0 * n5 + 100320856.0 / 1915538625.0 * n6;
        p.geod2auth[4] = 6007.0 / 14175.0 * n4 + 24496.0 / 467775.0 * n5 - 5884124.0 / 70945875.0 * n6;
        p.geod2auth[5] = -23356.0 / 66825.0 * n5 - 839792.0 / 19348875.0 * n6;
        p.geod2auth[6] = 570284222.0 / 1915538625.0 * n6;

        // AUTH2GEOD
        p.auth2geod[1] = 4.0 / 3.0 * n + 4.0 / 45.0 * n2 - 16.0 / 35.0 * n3 - 2582.0 / 14175.0 * n4 + 60136.0 / 467775.0 * n5 + 28112932.0 / 212837625.0 * n6;
        p.auth2geod[2] = 46.0 / 45.0 * n2 + 152.0 / 945.0 * n3 - 11966.0 / 14175.0 * n4 - 21016.0 / 51975.0 * n5 + 251310128.0 / 638512875.0 * n6;
        p.auth2geod[3] = 3044.0 / 2835.0 * n3 + 3802.0 / 14175.0 * n4 - 94388.0 / 66825.0 * n5 - 8797648.0 / 10945935.0 * n6;
        p.auth2geod[4] = 6059.0 / 4725.0 * n4 + 41072.0 / 93555.0 * n5 - 1472637812.0 / 638512875.0 * n6;
        p.auth2geod[5] = 768272.0 / 467775.0 * n5 + 455935736.0 / 638512875.0 * n6;
        p.auth2geod[6] = 4210684958.0 / 1915538625.0 * n6;

        // GEOD2GEOC
        p.geod2geoc[1] = -2.0 * n + 2.0 * n3 - 2.0 * n5;
        p.geod2geoc[2] = 2.0 * n2 - 4.0 * n4 + 6.0 * n6;
        p.geod2geoc[3] = -8.0 / 3.0 * n3 + 8.0 * n5;
        p.geod2geoc[4] = 4.0 * n4 - 16.0 * n6;
        p.geod2geoc[5] = -32.0 / 5.0 * n5;
        p.geod2geoc[6] = 32.0 / 3.0 * n6;

        // GEOC2GEOD
        p.geoc2geod[1] = 2.0 * n - 2.0 * n3 + 2.0 * n5;
        p.geoc2geod[2] = 2.0 * n2 - 4.0 * n4 + 6.0 * n6;
        p.geoc2geod[3] = 8.0 / 3.0 * n3 - 8.0 * n5;
        p.geoc2geod[4] = 4.0 * n4 - 16.0 * n6;
        p.geoc2geod[5] = 32.0 / 5.0 * n5;
        p.geoc2geod[6] = 32.0 / 3.0 * n6;

        // GEOD2PARA
        p.geod2para[1] = -n;
        p.geod2para[2] = 1.0 / 2.0 * n2;
        p.geod2para[3] = -1.0 / 3.0 * n3;
        p.geod2para[4] = 1.0 / 4.0 * n4;
        p.geod2para[5] = -1.0 / 5.0 * n5;
        p.geod2para[6] = 1.0 / 6.0 * n6;

        // PARA2GEOD
        p.para2geod[1] = n;
        p.para2geod[2] = 1.0 / 2.0 * n2;
        p.para2geod[3] = 1.0 / 3.0 * n3;
        p.para2geod[4] = 1.0 / 4.0 * n4;
        p.para2geod[5] = 1.0 / 5.0 * n5;
        p.para2geod[6] = 1.0 / 6.0 * n6;

        // CONF2RECT
        p.conf2rect[1] = 1.0 / 2.0 * n - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3 + 41.0 / 180.0 * n4 - 127.0 / 288.0 * n5 + 7891.0 / 37800.0 * n6;
        p.conf2rect[2] = 13.0 / 48.0 * n2 - 3.0 / 5.0 * n3 + 557.0 / 1440.0 * n4 + 281.0 / 630.0 * n5 - 1983433.0 / 1935360.0 * n6;
        p.conf2rect[3] = 61.0 / 240.0 * n3 - 103.0 / 140.0 * n4 + 15061.0 / 26880.0 * n5 + 167603.0 / 181440.0 * n6;
        p.conf2rect[4] = 49561.0 / 161280.0 * n4 - 179.0 / 168.0 * n5 + 6601661.0 / 7257600.0 * n6;
        p.conf2rect[5] = 34729.0 / 80640.0 * n5 - 3418889.0 / 1995840.0 * n6;
        p.conf2rect[6] = 212378941.0 / 319334400.0 * n6;

        // RECT2CONF
        p.rect2conf[1] = -1.0 / 2.0 * n + 2.0 / 3.0 * n2 - 37.0 / 96.0 * n3 + 1.0 / 360.0 * n4 + 81.0 / 512.0 * n5 - 96199.0 / 604800.0 * n6;
        p.rect2conf[2] = -1.0 / 48.0 * n2 - 1.0 / 15.0 * n3 + 437.0 / 1440.0 * n4 - 46.0 / 105.0 * n5 + 1118711.0 / 3870720.0 * n6;
        p.rect2conf[3] = -17.0 / 480.0 * n3 + 37.0 / 840.0 * n4 + 209.0 / 4480.0 * n5 - 5569.0 / 90720.0 * n6;
        p.rect2conf[4] = -4397.0 / 161280.0 * n4 + 11.0 / 504.0 * n5 + 830251.0 / 7257600.0 * n6;
        p.rect2conf[5] = -4583.0 / 161280.0 * n5 + 108847.0 / 3991680.0 * n6;
        p.rect2conf[6] = -20648693.0 / 638668800.0 * n6;

        // AUTH2RECT
        p.auth2rect[1] = -1.0 / 6.0 * n + 4.0 / 45.0 * n2 + 121.0 / 1680.0 * n3 - 1609.0 / 28350.0 * n4 - 384229.0 / 14968800.0 * n5 + 12674323.0 / 851350500.0 * n6;
        p.auth2rect[2] = -29.0 / 720.0 * n2 + 26.0 / 945.0 * n3 + 16463.0 / 453600.0 * n4 - 431.0 / 17325.0 * n5 - 31621753811.0 / 1307674368000.0 * n6;
        p.auth2rect[3] = -1003.0 / 45360.0 * n3 + 449.0 / 28350.0 * n4 + 3746047.0 / 119750400.0 * n5 - 32844781.0 / 1751349600.0 * n6;
        p.auth2rect[4] = -40457.0 / 2419200.0 * n4 + 629.0 / 53460.0 * n5 + 10650637121.0 / 326918592000.0 * n6;
        p.auth2rect[5] = -1800439.0 / 119750400.0 * n5 + 205072597.0 / 20432412000.0 * n6;
        p.auth2rect[6] = -59109051671.0 / 3923023104000.0 * n6;

        // RECT2AUTH
        p.rect2auth[1] = 1.0 / 6.0 * n - 4.0 / 45.0 * n2 - 817.0 / 10080.0 * n3 + 1297.0 / 18900.0 * n4 + 7764059.0 / 239500800.0 * n5 - 9292991.0 / 302702400.0 * n6;
        p.rect2auth[2] = 49.0 / 720.0 * n2 - 2.0 / 35.0 * n3 - 29609.0 / 453600.0 * n4 + 35474.0 / 467775.0 * n5 + 36019108271.0 / 871782912000.0 * n6;
        p.rect2auth[3] = 4463.0 / 90720.0 * n3 - 2917.0 / 56700.0 * n4 - 4306823.0 / 59875200.0 * n5 + 3026004511.0 / 30648618000.0 * n6;
        p.rect2auth[4] = 331799.0 / 7257600.0 * n4 - 102293.0 / 1871100.0 * n5 - 368661577.0 / 4036032000.0 * n6;
        p.rect2auth[5] = 11744233.0 / 239500800.0 * n5 - 875457073.0 / 13621608000.0 * n6;
        p.rect2auth[6] = 453002260127.0 / 7846046208000.0 * n6;

        // GEOC2RECT
        p.geoc2rect[1] = 1.0 / 2.0 * n + 13.0 / 16.0 * n3 - 15.0 / 32.0 * n5;
        p.geoc2rect[2] = -1.0 / 16.0 * n2 + 33.0 / 32.0 * n4 - 1673.0 / 2048.0 * n6;
        p.geoc2rect[3] = -5.0 / 16.0 * n3 + 349.0 / 256.0 * n5;
        p.geoc2rect[4] = -261.0 / 512.0 * n4 + 963.0 / 512.0 * n6;
        p.geoc2rect[5] = -921.0 / 1280.0 * n5;
        p.geoc2rect[6] = -6037.0 / 6144.0 * n6;

        // RECT2GEOC
        p.rect2geoc[1] = -1.0 / 2.0 * n - 23.0 / 32.0 * n3 + 499.0 / 1536.0 * n5;
        p.rect2geoc[2] = 5.0 / 16.0 * n2 - 5.0 / 96.0 * n4 + 6565.0 / 12288.0 * n6;
        p.rect2geoc[3] = 1.0 / 32.0 * n3 - 77.0 / 128.0 * n5;
        p.rect2geoc[4] = 283.0 / 1536.0 * n4 - 4037.0 / 7680.0 * n6;
        p.rect2geoc[5] = 1301.0 / 7680.0 * n5;
        p.rect2geoc[6] = 17089.0 / 61440.0 * n6;

        // PARA2RECT
        p.para2rect[1] = -1.0 / 2.0 * n + 3.0 / 16.0 * n3 - 1.0 / 32.0 * n5;
        p.para2rect[2] = -1.0 / 16.0 * n2 + 1.0 / 32.0 * n4 - 9.0 / 2048.0 * n6;
        p.para2rect[3] = -1.0 / 48.0 * n3 + 3.0 / 256.0 * n5;
        p.para2rect[4] = -5.0 / 512.0 * n4 + 3.0 / 512.0 * n6;
        p.para2rect[5] = -7.0 / 1280.0 * n5;
        p.para2rect[6] = -7.0 / 2048.0 * n6;

        // RECT2PARA
        p.rect2para[1] = 1.0 / 2.0 * n - 9.0 / 32.0 * n3 + 205.0 / 1536.0 * n5;
        p.rect2para[2] = 5.0 / 16.0 * n2 - 37.0 / 96.0 * n4 + 1335.0 / 4096.0 * n6;
        p.rect2para[3] = 29.0 / 96.0 * n3 - 75.0 / 128.0 * n5;
        p.rect2para[4] = 539.0 / 1536.0 * n4 - 2391.0 / 2560.0 * n6;
        p.rect2para[5] = 3467.0 / 7680.0 * n5;
        p.rect2para[6] = 38081.0 / 61440.0 * n6;

        // AUTH2CONF
        p.auth2conf[1] = -2.0 / 3.0 * n + 34.0 / 45.0 * n2 - 88.0 / 315.0 * n3 - 2312.0 / 14175.0 * n4 + 27128.0 / 93555.0 * n5 - 55271278.0 / 212837625.0 * n6;
        p.auth2conf[2] = 1.0 / 45.0 * n2 - 184.0 / 945.0 * n3 + 6079.0 / 14175.0 * n4 - 65864.0 / 155925.0 * n5 + 106691108.0 / 638512875.0 * n6;
        p.auth2conf[3] = -106.0 / 2835.0 * n3 + 772.0 / 14175.0 * n4 - 14246.0 / 467775.0 * n5 + 5921152.0 / 54729675.0 * n6;
        p.auth2conf[4] = -167.0 / 9450.0 * n4 - 5312.0 / 467775.0 * n5 + 75594328.0 / 638512875.0 * n6;
        p.auth2conf[5] = -248.0 / 13365.0 * n5 + 2837636.0 / 638512875.0 * n6;
        p.auth2conf[6] = -34761247.0 / 1915538625.0 * n6;

        // CONF2AUTH
        p.conf2auth[1] = 2.0 / 3.0 * n - 34.0 / 45.0 * n2 + 46.0 / 315.0 * n3 + 2458.0 / 4725.0 * n4 - 55222.0 / 93555.0 * n5 + 2706758.0 / 42567525.0 * n6;
        p.conf2auth[2] = 19.0 / 45.0 * n2 - 256.0 / 315.0 * n3 + 3413.0 / 14175.0 * n4 + 516944.0 / 467775.0 * n5 - 340492279.0 / 212837625.0 * n6;
        p.conf2auth[3] = 248.0 / 567.0 * n3 - 15958.0 / 14175.0 * n4 + 206834.0 / 467775.0 * n5 + 4430783356.0 / 1915538625.0 * n6;
        p.conf2auth[4] = 16049.0 / 28350.0 * n4 - 832976.0 / 467775.0 * n5 + 62016436.0 / 70945875.0 * n6;
        p.conf2auth[5] = 15602.0 / 18711.0 * n5 - 651151712.0 / 212837625.0 * n6;
        p.conf2auth[6] = 2561772812.0 / 1915538625.0 * n6;

        // GEOC2CONF
        p.geoc2conf[1] = 2.0 / 3.0 * n2 + 2.0 / 3.0 * n3 - 2.0 / 9.0 * n4 - 14.0 / 45.0 * n5 + 1042.0 / 4725.0 * n6;
        p.geoc2conf[2] = -1.0 / 3.0 * n2 + 4.0 / 15.0 * n3 + 43.0 / 45.0 * n4 - 4.0 / 45.0 * n5 - 712.0 / 945.0 * n6;
        p.geoc2conf[3] = -2.0 / 5.0 * n3 + 2.0 / 105.0 * n4 + 124.0 / 105.0 * n5 + 274.0 / 2835.0 * n6;
        p.geoc2conf[4] = -55.0 / 126.0 * n4 - 16.0 / 105.0 * n5 + 21068.0 / 14175.0 * n6;
        p.geoc2conf[5] = -22.0 / 45.0 * n5 - 9202.0 / 31185.0 * n6;
        p.geoc2conf[6] = -90263.0 / 155925.0 * n6;

        // CONF2GEOC
        p.conf2geoc[1] = -2.0 / 3.0 * n2 - 2.0 / 3.0 * n3 + 4.0 / 9.0 * n4 + 2.0 / 9.0 * n5 - 3658.0 / 4725.0 * n6;
        p.conf2geoc[2] = 1.0 / 3.0 * n2 - 4.0 / 15.0 * n3 - 23.0 / 45.0 * n4 + 68.0 / 45.0 * n5 + 61.0 / 135.0 * n6;
        p.conf2geoc[3] = 2.0 / 5.0 * n3 - 24.0 / 35.0 * n4 - 46.0 / 35.0 * n5 + 9446.0 / 2835.0 * n6;
        p.conf2geoc[4] = 83.0 / 126.0 * n4 - 80.0 / 63.0 * n5 - 34712.0 / 14175.0 * n6;
        p.conf2geoc[5] = 52.0 / 45.0 * n5 - 2362.0 / 891.0 * n6;
        p.conf2geoc[6] = 335882.0 / 155925.0 * n6;

        // PARA2CONF
        p.para2conf[1] = -n + 2.0 / 3.0 * n2 - 16.0 / 45.0 * n4 + 2.0 / 5.0 * n5 - 998.0 / 4725.0 * n6;
        p.para2conf[2] = 1.0 / 6.0 * n2 - 2.0 / 5.0 * n3 + 19.0 / 45.0 * n4 - 22.0 / 105.0 * n5 - 2.0 / 27.0 * n6;
        p.para2conf[3] = -1.0 / 15.0 * n3 + 16.0 / 105.0 * n4 - 22.0 / 105.0 * n5 + 116.0 / 567.0 * n6;
        p.para2conf[4] = 17.0 / 1260.0 * n4 - 8.0 / 105.0 * n5 + 2123.0 / 14175.0 * n6;
        p.para2conf[5] = -1.0 / 105.0 * n5 + 128.0 / 4455.0 * n6;
        p.para2conf[6] = 149.0 / 311850.0 * n6;

        // CONF2PARA
        p.conf2para[1] = n - 2.0 / 3.0 * n2 - 1.0 / 3.0 * n3 + 38.0 / 45.0 * n4 - 1.0 / 3.0 * n5 - 3118.0 / 4725.0 * n6;
        p.conf2para[2] = 5.0 / 6.0 * n2 - 14.0 / 15.0 * n3 - 7.0 / 9.0 * n4 + 50.0 / 21.0 * n5 - 247.0 / 270.0 * n6;
        p.conf2para[3] = 16.0 / 15.0 * n3 - 34.0 / 21.0 * n4 - 5.0 / 3.0 * n5 + 17564.0 / 2835.0 * n6;
        p.conf2para[4] = 2069.0 / 1260.0 * n4 - 28.0 / 9.0 * n5 - 49877.0 / 14175.0 * n6;
        p.conf2para[5] = 883.0 / 315.0 * n5 - 28244.0 / 4455.0 * n6;
        p.conf2para[6] = 797222.0 / 155925.0 * n6;

        // GEOC2AUTH
        p.geoc2auth[1] = 2.0 / 3.0 * n - 4.0 / 45.0 * n2 + 62.0 / 105.0 * n3 + 778.0 / 4725.0 * n4 - 193082.0 / 467775.0 * n5 - 4286228.0 / 42567525.0 * n6;
        p.geoc2auth[2] = 4.0 / 45.0 * n2 - 32.0 / 315.0 * n3 + 12338.0 / 14175.0 * n4 + 92696.0 / 467775.0 * n5 - 61623938.0 / 70945875.0 * n6;
        p.geoc2auth[3] = -524.0 / 2835.0 * n3 - 1618.0 / 14175.0 * n4 + 612536.0 / 467775.0 * n5 + 427003576.0 / 1915538625.0 * n6;
        p.geoc2auth[4] = -5933.0 / 14175.0 * n4 - 8324.0 / 66825.0 * n5 + 427770788.0 / 212837625.0 * n6;
        p.geoc2auth[5] = -320044.0 / 467775.0 * n5 - 9153184.0 / 70945875.0 * n6;
        p.geoc2auth[6] = -1978771378.0 / 1915538625.0 * n6;

        // AUTH2GEOC
        p.auth2geoc[1] = -2.0 / 3.0 * n + 4.0 / 45.0 * n2 - 158.0 / 315.0 * n3 - 2102.0 / 14175.0 * n4 + 109042.0 / 467775.0 * n5 + 216932.0 / 2627625.0 * n6;
        p.auth2geoc[2] = 16.0 / 45.0 * n2 - 16.0 / 945.0 * n3 + 934.0 / 14175.0 * n4 - 7256.0 / 155925.0 * n5 + 117952358.0 / 638512875.0 * n6;
        p.auth2geoc[3] = -232.0 / 2835.0 * n3 + 922.0 / 14175.0 * n4 - 25286.0 / 66825.0 * n5 - 7391576.0 / 54729675.0 * n6;
        p.auth2geoc[4] = 719.0 / 4725.0 * n4 + 268.0 / 18711.0 * n5 - 67048172.0 / 638512875.0 * n6;
        p.auth2geoc[5] = 14354.0 / 467775.0 * n5 + 46774256.0 / 638512875.0 * n6;
        p.auth2geoc[6] = 253129538.0 / 1915538625.0 * n6;

        // PARA2AUTH
        p.para2auth[1] = -1.0 / 3.0 * n - 4.0 / 45.0 * n2 + 32.0 / 315.0 * n3 + 34.0 / 675.0 * n4 + 2476.0 / 467775.0 * n5 - 70496.0 / 8513505.0 * n6;
        p.para2auth[2] = -7.0 / 90.0 * n2 - 4.0 / 315.0 * n3 + 74.0 / 2025.0 * n4 + 3992.0 / 467775.0 * n5 + 53836.0 / 212837625.0 * n6;
        p.para2auth[3] = -83.0 / 2835.0 * n3 + 2.0 / 14175.0 * n4 + 7052.0 / 467775.0 * n5 - 661844.0 / 1915538625.0 * n6;
        p.para2auth[4] = -797.0 / 56700.0 * n4 + 934.0 / 467775.0 * n5 + 1425778.0 / 212837625.0 * n6;
        p.para2auth[5] = -3673.0 / 467775.0 * n5 + 390088.0 / 212837625.0 * n6;
        p.para2auth[6] = -18623681.0 / 3831077250.0 * n6;

        // AUTH2PARA
        p.auth2para[1] = 1.0 / 3.0 * n + 4.0 / 45.0 * n2 - 46.0 / 315.0 * n3 - 1082.0 / 14175.0 * n4 + 11824.0 / 467775.0 * n5 + 7947332.0 / 212837625.0 * n6;
        p.auth2para[2] = 17.0 / 90.0 * n2 + 68.0 / 945.0 * n3 - 338.0 / 2025.0 * n4 - 16672.0 / 155925.0 * n5 + 39946703.0 / 638512875.0 * n6;
        p.auth2para[3] = 461.0 / 2835.0 * n3 + 1102.0 / 14175.0 * n4 - 101069.0 / 467775.0 * n5 - 255454.0 / 1563705.0 * n6;
        p.auth2para[4] = 3161.0 / 18900.0 * n4 + 1786.0 / 18711.0 * n5 - 189032762.0 / 638512875.0 * n6;
        p.auth2para[5] = 88868.0 / 467775.0 * n5 + 80274086.0 / 638512875.0 * n6;
        p.auth2para[6] = 880980241.0 / 3831077250.0 * n6;

        // PARA2GEOC
        p.para2geoc[1] = -n;
        p.para2geoc[2] = 1.0 / 2.0 * n2;
        p.para2geoc[3] = -1.0 / 3.0 * n3;
        p.para2geoc[4] = 1.0 / 4.0 * n4;
        p.para2geoc[5] = -1.0 / 5.0 * n5;
        p.para2geoc[6] = 1.0 / 6.0 * n6;

        // GEOC2PARA
        p.geoc2para[1] = n;
        p.geoc2para[2] = 1.0 / 2.0 * n2;
        p.geoc2para[3] = 1.0 / 3.0 * n3;
        p.geoc2para[4] = 1.0 / 4.0 * n4;
        p.geoc2para[5] = 1.0 / 5.0 * n5;
        p.geoc2para[6] = 1.0 / 6.0 * n6;

        p
    }
}

// ============================================================
// FUNCIONES DE CONVERSIÓN DE LATITUD
fn lat2latn6(a: &[f64; 7], lat1: f64) -> f64 {
    let cp = (2.0 * lat1).cos();
    lat1 + (cp * (cp * (cp * (cp * (32.0 * a[6] * cp + 16.0 * a[5]) - 32.0 * a[6] + 8.0 * a[4]) - 12.0 * a[5] + 4.0 * a[3]) + 6.0 * a[6] - 4.0 * a[4] + 2.0 * a[2]) + a[5] - a[3] + a[1]) * (2.0 * lat1).sin()
}

fn partial_from_latn6(a: &[f64; 7], phi: f64) -> f64 {
    let cp = phi.cos();
    let cp2 = cp * cp;
    let latp = cp2
        * (cp2
            * (cp2 * (cp2 * (cp2 * (24576.0 * a[6] * cp2 - 73728.0 * a[6] + 5120.0 * a[5]) + 82944.0 * a[6] - 12800.0 * a[5] + 1024.0 * a[4]) - 43008.0 * a[6] + 11200.0 * a[5] - 2048.0 * a[4] + 192.0 * a[3])
                + 10080.0 * a[6]
                - 4000.0 * a[5]
                + 1280.0 * a[4]
                - 288.0 * a[3]
                + 32.0 * a[2])
            - 864.0 * a[6]
            + 500.0 * a[5]
            - 256.0 * a[4]
            + 108.0 * a[3]
            - 32.0 * a[2]
            + 4.0 * a[1])
        + 12.0 * a[6]
        - 10.0 * a[5]
        + 8.0 * a[4]
        - 6.0 * a[3]
        + 4.0 * a[2]
        - 2.0 * a[1];
    1.0 + latp
}

// ============================================================
// FUNCIONES GEODÉSICAS BÁSICAS
fn radios_ell(c: f64, ep: f64, latitude: f64) -> (f64, f64, f64, f64) {
    let cosphi = latitude.cos();
    let v = (1.0 + (ep * cosphi).powi(2)).sqrt();
    let rn = c / v;
    let rm = c / (v * v * v);
    let rg = c / (v * v);
    let r = rn * cosphi;
    (rn, rm, rg, r)
}

fn phi2psi(elli: &Pelipsoide, phi: f64) -> f64 {
    let (sphi, cphi) = sincos(phi);
    let u = sphi * (1.0 - elli.e2);
    toquadcirc(arg(cphi, u))
}

fn psi2phi(elli: &Pelipsoide, psi: f64) -> f64 {
    let (spsi, cpsi) = sincos(psi);
    let v = cpsi * (1.0 - elli.e2);
    toquadcirc(arg(v, spsi))
}

fn phi2beta(elli: &Pelipsoide, phi: f64) -> f64 {
    let (sphi, cphi) = sincos(phi);
    toquadcirc(arg(cphi, (1.0 - elli.f) * sphi))
}

fn beta2phi(elli: &Pelipsoide, beta: f64) -> f64 {
    let (sbeta, cbeta) = sincos(beta);
    toquadcirc(arg((1.0 - elli.f) * cbeta, sbeta))
}

fn beta_partial_phi(elli: &Pelipsoide, phi: f64) -> f64 {
    (1.0 - elli.f) / (1.0 - elli.f * (2.0 - elli.f) * phi.sin().powi(2))
}

// ============================================================
// CURVA DE ALINEACIÓN (ALIGN)
fn cs_align(elli: &Pelipsoide, beta1: f64, l1: f64, beta2: f64, l2: f64) -> [f64; 6] {
    let (sbeta1, cbeta1) = sincos(beta1);
    let (sl1, cl1) = sincos(l1);
    let (sbeta2, cbeta2) = sincos(beta2);
    let (sl2, cl2) = sincos(l2);

    let x1 = elli.a * cbeta1 * cl1;
    let y1 = elli.a * cbeta1 * sl1;
    let z1 = elli.b * sbeta1;
    let x2 = elli.a * cbeta2 * cl2;
    let y2 = elli.a * cbeta2 * sl2;
    let z2 = elli.b * sbeta2;

    let ab = elli.a / elli.b;
    let ba = elli.b / elli.a;

    let c1 = y1 * ab - y1 * ba - y2 * ab + y2 * ba;
    let c2 = -x1 * ab + x1 * ba + x2 * ab - x2 * ba;
    let c3 = (y1 / elli.a) * z2 - (y2 / elli.a) * z1;
    let c4 = (-x1 / elli.a) * z2 + (x2 / elli.a) * z1;
    let c5 = (x1 / elli.b) * y2 - (x2 / elli.b) * y1;

    [0.0, c1, c2, c3, c4, c5]
}

fn g_implicita_align(cs: &[f64; 6], beta: f64, longitude: f64) -> f64 {
    let (sbeta, cbeta) = sincos(beta);
    let (sl, cl) = sincos(longitude);
    cs[1] * cbeta * sbeta * cl + cs[2] * cbeta * sbeta * sl + cs[3] * cbeta * cl + cs[4] * cbeta * sl + cs[5] * sbeta
}

fn g_partial_align(cs: &[f64; 6], beta: f64, longitude: f64) -> (f64, f64) {
    let c2beta = (2.0 * beta).cos();
    let (sbeta, cbeta) = sincos(beta);
    let (sl, cl) = sincos(longitude);

    let g_beta = c2beta * (cs[1] * cl + cs[2] * sl) + sbeta * (-cs[3] * cl - cs[4] * sl) + cs[5] * cbeta;
    let g_l = cbeta * (sbeta * (cs[2] * cl - cs[1] * sl) - cs[3] * sl + cs[4] * cl);

    (g_beta, g_l)
}

fn acimut_align(elli: &Pelipsoide, cs: &[f64; 6], beta: f64, longitude: f64) -> f64 {
    let cbeta = beta.cos();
    let (g_beta, g_l) = g_partial_align(cs, beta, longitude);
    let y_tanalpha = cbeta * g_beta;
    let x_tanalpha = -(1.0 - elli.e2 * cbeta * cbeta).sqrt() * g_l;
    arg(x_tanalpha, y_tanalpha)
}

fn l2beta_align(cs: &[f64; 6], mut beta: f64, l: f64) -> f64 {
    let mut sl: f64;
    let mut cl: f64;
    (sl, cl) = sincos(l);

    for _ in 0..30 {
        let c2beta = (2.0 * beta).cos();
        let (mut sbeta, cbeta) = sincos(beta);

        let g_beta = c2beta * (cs[1] * cl + cs[2] * sl) + sbeta * (-cs[3] * cl - cs[4] * sl) + cs[5] * cbeta;
        //let g_l = cbeta * (sbeta * (cs[2] * cl - cs[1] * sl) - cs[3] * sl + cs[4] * cl);
        let g = cs[1] * cbeta * sbeta * cl + cs[2] * cbeta * sbeta * sl + cs[3] * cbeta * cl + cs[4] * cbeta * sl + cs[5] * sbeta;

        let mut delta_beta = g / (g_beta + sgn(g_beta) * EPSILON_MAQ);
        if delta_beta.abs() > MAXSTEP_ANG {
            delta_beta = MAXSTEP_ANG.copysign(delta_beta);
        }
        beta -= delta_beta;

        if delta_beta.abs() < EPSILON_ANG {
            break;
        }
        sbeta = beta.sin();
        beta = sbeta.asin();
        (sl, cl) = sincos(l);
    }
    toquadcirc(beta)
}

fn beta2l_align(cs: &[f64; 6], beta: f64, mut l: f64) -> f64 {
    let (mut sl, mut cl) = sincos(l);

    for _ in 0..30 {
        //let c2beta = (2.0 * beta).cos();
        let (sbeta, cbeta) = sincos(beta);

        //let g_beta = c2beta * (cs[1] * cl + cs[2] * sl) + sbeta * (-cs[3] * cl - cs[4] * sl) + cs[5] * cbeta;
        let g_l = cbeta * (sbeta * (cs[2] * cl - cs[1] * sl) - cs[3] * sl + cs[4] * cl);
        let g = cs[1] * cbeta * sbeta * cl + cs[2] * cbeta * sbeta * sl + cs[3] * cbeta * cl + cs[4] * cbeta * sl + cs[5] * sbeta;

        let mut delta_l = g / (g_l + sgn(g_l) * EPSILON_MAQ);
        if delta_l.abs() > MAXSTEP_ANG {
            delta_l = MAXSTEP_ANG.copysign(delta_l);
        }
        l -= delta_l;
        (sl, cl) = sincos(l);

        if delta_l.abs() < EPSILON_ANG {
            break;
        }
    }
    tosemicirc(l)
}

fn calc_beta0_align(cs: &[f64; 6], mut beta: f64, mut l: f64) -> (f64, f64) {
    let (mut sbeta, mut cbeta) = sincos(beta);
    let (mut sl, mut cl) = sincos(l);

    for _ in 0..30 {
        let c2beta = (2.0 * beta).cos();
        let g = cs[1] * cbeta * sbeta * cl + cs[2] * cbeta * sbeta * sl + cs[3] * cbeta * cl + cs[4] * cbeta * sl + cs[5] * sbeta;
        let g_beta = c2beta * (cs[1] * cl + cs[2] * sl) + sbeta * (-cs[3] * cl - cs[4] * sl) + cs[5] * cbeta;
        let g_l = cbeta * (sbeta * (cs[2] * cl - cs[1] * sl) - cs[3] * sl + cs[4] * cl);
        let g_l_b = 2.0 * cbeta * cbeta * (cs[2] * cl - cs[1] * sl) + sbeta * (cs[3] * sl - cs[4] * cl) + cs[1] * sl - cs[2] * cl;
        let g_l_l = -cbeta * (sbeta * (cs[1] * cl + cs[2] * sl) + cs[3] * cl + cs[4] * sl);

        let a = g_beta;
        let b = g_l;
        let c = g_l_b;
        let d = g_l_l;
        let f1 = g;
        let f2 = g_l;
        let det = a * d - b * c;
        let sgndet = sgn(det);

        let mut delta_beta = (d * f1 - b * f2) / (det + sgndet * EPSILON_MAQ);
        let mut delta_l = (-c * f1 + a * f2) / (det + sgndet * EPSILON_MAQ);

        if delta_l.abs() > MAXSTEP_ANG {
            delta_l = MAXSTEP_ANG.copysign(delta_l);
        }
        if delta_beta.abs() > MAXSTEP_ANG {
            delta_beta = MAXSTEP_ANG.copysign(delta_beta);
        }

        beta -= delta_beta;
        l -= delta_l;
        (sbeta, cbeta) = sincos(beta);
        (sl, cl) = sincos(l);

        if delta_beta.abs().max(delta_l.abs()) < EPSILON_ANG {
            break;
        }
    }

    beta = sbeta.asin();
    l = arg(cl, sl);
    (beta, l)
}

// ============================================================
// SECCIÓN CENTRAL (CENTRAL)
fn cs_central(elli: &Pelipsoide, beta1: f64, l1: f64, beta2: f64, l2: f64) -> [f64; 4] {
    let (sbeta1, cbeta1) = sincos(beta1);
    let (sl1, cl1) = sincos(l1);
    let (sbeta2, cbeta2) = sincos(beta2);
    let (sl2, cl2) = sincos(l2);

    let x1 = elli.a * cbeta1 * cl1;
    let y1 = elli.a * cbeta1 * sl1;
    let z1 = elli.b * sbeta1;
    let x2 = elli.a * cbeta2 * cl2;
    let y2 = elli.a * cbeta2 * sl2;
    let z2 = elli.b * sbeta2;

    let c1 = (y1 / elli.b) * z2 - y2 * (z1 / elli.b);
    let c2 = -(x1 / elli.b) * z2 + x2 * (z1 / elli.b);
    let c3 = (x1 / elli.a) * y2 - x2 * (y1 / elli.a);

    [0.0, c1, c2, c3]
}

fn g_implicita_central(cs: &[f64; 4], beta: f64, longitude: f64) -> f64 {
    let (sbeta, cbeta) = sincos(beta);
    let (sl, cl) = sincos(longitude);
    cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta
}

fn g_partial_central(cs: &[f64; 4], beta: f64, longitude: f64) -> (f64, f64) {
    let (sbeta, cbeta) = sincos(beta);
    let (sl, cl) = sincos(longitude);

    let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
    let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;

    (g_beta, g_l)
}

fn acimut_central(elli: &Pelipsoide, cs: &[f64; 4], beta: f64, longitude: f64) -> f64 {
    let cbeta = beta.cos();
    let (g_beta, g_l) = g_partial_central(cs, beta, longitude);
    let y_tanalpha = cbeta * g_beta;
    let x_tanalpha = -(1.0 - elli.e2 * cbeta * cbeta).sqrt() * g_l;
    arg(x_tanalpha, y_tanalpha)
}

fn l2beta_central(cs: &[f64; 4], mut beta: f64, l: f64) -> f64 {
    let (sl, cl) = sincos(l);

    for _ in 0..30 {
        let (sbeta, cbeta) = sincos(beta);
        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        //let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let g = cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta;

        let mut delta_beta = g / (g_beta + sgn(g_beta) * EPSILON_MAQ);
        if delta_beta.abs() > MAXSTEP_ANG {
            delta_beta = MAXSTEP_ANG.copysign(delta_beta);
        }
        beta -= delta_beta;

        if delta_beta.abs() < EPSILON_ANG {
            break;
        }
        beta = toquadcirc(beta);
    }
    beta
}

fn beta2l_central(cs: &[f64; 4], beta: f64, mut l: f64) -> f64 {
    let (mut sl, mut cl) = sincos(l);

    for _ in 0..30 {
        let (sbeta, cbeta) = sincos(beta);
        //let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let g = cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta;

        let mut delta_l = g / (g_l + sgn(g_l) * EPSILON_MAQ);
        if delta_l.abs() > MAXSTEP_ANG {
            delta_l = MAXSTEP_ANG.copysign(delta_l);
        }
        l -= delta_l;
        (sl, cl) = sincos(l);

        if delta_l.abs() < EPSILON_ANG {
            break;
        }
    }
    tosemicirc(l)
}

fn calc_beta0_central(cs: &[f64; 4], mut beta: f64, mut l: f64) -> (f64, f64) {
    let (mut sbeta, mut cbeta) = sincos(beta);
    let (mut sl, mut cl) = sincos(l);

    for _ in 0..30 {
        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let g = cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta;
        let g_l_b = cs[1] * sbeta * sl - cs[2] * sbeta * cl;
        let g_l_l = -cs[1] * cbeta * cl - cs[2] * cbeta * sl;

        let a = g_beta;
        let b = g_l;
        let c = g_l_b;
        let d = g_l_l;
        let f1 = g;
        let f2 = g_l;
        let det = a * d - b * c;
        let sgndet = sgn(det);

        let mut delta_beta = (d * f1 - b * f2) / (det + sgndet * EPSILON_MAQ);
        let mut delta_l = (-c * f1 + a * f2) / (det + sgndet * EPSILON_MAQ);

        if delta_l.abs() > MAXSTEP_ANG {
            delta_l = MAXSTEP_ANG.copysign(delta_l);
        }
        if delta_beta.abs() > MAXSTEP_ANG {
            delta_beta = MAXSTEP_ANG.copysign(delta_beta);
        }

        beta -= delta_beta;
        l -= delta_l;
        (sbeta, cbeta) = sincos(beta);
        (sl, cl) = sincos(l);

        if delta_beta.abs().max(delta_l.abs()) < EPSILON_ANG {
            break;
        }
    }

    beta = sbeta.asin();
    l = arg(cl, sl);
    (beta, l)
}

// ============================================================
// SECCIÓN NORMAL (NORMAL)
fn cs_normal(elli: &Pelipsoide, beta1: f64, l1: f64, beta2: f64, l2: f64) -> [f64; 5] {
    let (sbeta1, cbeta1) = sincos(beta1);
    let (sl1, cl1) = sincos(l1);
    let (sbeta2, cbeta2) = sincos(beta2);
    let (sl2, cl2) = sincos(l2);

    let x1 = elli.a * cbeta1 * cl1;
    let y1 = elli.a * cbeta1 * sl1;
    let z1 = elli.b * sbeta1;
    let x2 = elli.a * cbeta2 * cl2;
    let y2 = elli.a * cbeta2 * sl2;
    let z2 = elli.b * sbeta2;

    let c1 = elli.a * (y1 / elli.b) * (z1 / elli.b) - (y1 / elli.a) * z1 + (y1 / elli.a) * z2 - elli.a * (y2 / elli.b) * (z1 / elli.b);
    let c2 = -elli.a * (x1 / elli.b) * (z1 / elli.b) + (x1 / elli.a) * z1 - (x1 / elli.a) * z2 + elli.a * (x2 / elli.b) * (z1 / elli.b);
    let c3 = elli.b * (x1 / elli.a) * (y2 / elli.a) - elli.b * (x2 / elli.a) * (y1 / elli.a);
    let c4 = (x1 / elli.b) * (y2 / elli.b) * z1 - x1 * (y2 / elli.a) * (z1 / elli.a) - (x2 / elli.b) * y1 * (z1 / elli.b) + x2 * (y1 / elli.a) * (z1 / elli.a);

    [0.0, c1, c2, c3, c4]
}

fn g_implicita_normal(cs: &[f64; 5], beta: f64, longitude: f64) -> f64 {
    let (sbeta, cbeta) = sincos(beta);
    let (sl, cl) = sincos(longitude);
    cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta + cs[4]
}

fn g_partial_normal(cs: &[f64; 5], beta: f64, longitude: f64) -> (f64, f64) {
    let (sbeta, cbeta) = sincos(beta);
    let (sl, cl) = sincos(longitude);

    let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
    let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;

    (g_beta, g_l)
}

fn acimut_normal(elli: &Pelipsoide, cs: &[f64; 5], beta: f64, longitude: f64) -> f64 {
    let cbeta = beta.cos();
    let (g_beta, g_l) = g_partial_normal(cs, beta, longitude);
    let y_tanalpha = cbeta * g_beta;
    let x_tanalpha = -(1.0 - elli.e2 * cbeta * cbeta).sqrt() * g_l;
    arg(x_tanalpha, y_tanalpha)
}

fn l2beta_normal(cs: &[f64; 5], mut beta: f64, l: f64) -> f64 {
    let (sl, cl) = sincos(l);

    for _ in 0..30 {
        let (sbeta, cbeta) = sincos(beta);
        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        //let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let g = cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta + cs[4];

        let mut delta_beta = g / (g_beta + sgn(g_beta) * EPSILON_MAQ);
        if delta_beta.abs() > MAXSTEP_ANG {
            delta_beta = MAXSTEP_ANG.copysign(delta_beta);
        }
        beta -= delta_beta;

        if delta_beta.abs() < EPSILON_ANG {
            break;
        }
        beta = toquadcirc(beta);
    }
    beta
}

fn beta2l_normal(cs: &[f64; 5], beta: f64, mut l: f64) -> f64 {
    let (mut sl, mut cl) = sincos(l);

    for _ in 0..30 {
        let (sbeta, cbeta) = sincos(beta);
        //let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let g = cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta + cs[4];

        let mut delta_l = g / (g_l + sgn(g_l) * EPSILON_MAQ);
        if delta_l.abs() > MAXSTEP_ANG {
            delta_l = MAXSTEP_ANG.copysign(delta_l);
        }
        l -= delta_l;
        (sl, cl) = sincos(l);

        if delta_l.abs() < EPSILON_ANG {
            break;
        }
    }
    l
}

fn calc_beta0_normal(cs: &[f64; 5], mut beta: f64, mut l: f64) -> (f64, f64) {
    let (mut sbeta, mut cbeta) = sincos(beta);
    let (mut sl, mut cl) = sincos(l);

    for _ in 0..30 {
        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let g = cs[1] * cbeta * cl + cs[2] * cbeta * sl + cs[3] * sbeta + cs[4];
        let g_l_b = cs[1] * sbeta * sl - cs[2] * sbeta * cl;
        let g_l_l = -cs[1] * cbeta * cl - cs[2] * cbeta * sl;

        let a = g_beta;
        let b = g_l;
        let c = g_l_b;
        let d = g_l_l;
        let f1 = g;
        let f2 = g_l;
        let det = a * d - b * c;
        let sgndet = sgn(det);

        let mut delta_beta = (d * f1 - b * f2) / (det + sgndet * EPSILON_MAQ);
        let mut delta_l = (-c * f1 + a * f2) / (det + sgndet * EPSILON_MAQ);

        if delta_l.abs() > MAXSTEP_ANG {
            delta_l = MAXSTEP_ANG.copysign(delta_l);
        }
        if delta_beta.abs() > MAXSTEP_ANG {
            delta_beta = MAXSTEP_ANG.copysign(delta_beta);
        }

        beta -= delta_beta;
        l -= delta_l;
        (sbeta, cbeta) = sincos(beta);
        (sl, cl) = sincos(l);

        if delta_beta.abs().max(delta_l.abs()) < EPSILON_ANG {
            break;
        }
    }

    beta = sbeta.asin();
    l = arg(cl, sl);
    (beta, l)
}

// ============================================================
// ECUACIONES DIFERENCIALES Y SOLVER RK45
const A_RK45: [f64; 7] = [0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const B_RK45: [[f64; 7]; 7] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0],
];
const C5_RK45: [f64; 7] = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0];
const C4_RK45: [f64; 7] = [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0];

fn rk45_step(y: &mut [f64], derivs: impl Fn(&[f64], f64) -> Vec<f64>, h: f64, x: f64) -> (Vec<f64>, f64, f64) {
    let n = y.len();
    let mut k = vec![vec![0.0; n]; 7];

    k[0] = derivs(y, x);
    for i in 1..7 {
        let xi = x + A_RK45[i] * h;
        let mut yi = y.to_vec();
        for j in 0..n {
            for m in 0..i {
                yi[j] += h * B_RK45[i][m] * k[m][j];
            }
        }
        k[i] = derivs(&yi, xi);
    }

    let mut y5 = y.to_vec();
    let mut y4 = y.to_vec();
    for j in 0..n {
        for i in 0..7 {
            y5[j] += h * C5_RK45[i] * k[i][j];
            y4[j] += h * C4_RK45[i] * k[i][j];
        }
    }

    let error0 = (y5[0] - y4[0]).abs();
    let error1 = (y5[1] - y4[1]).abs();
    (y5, error1, error0)
}

// ============================================================
// FUNCIONES PARA EL PROBLEMA DIRECTO (CentralSect2GEO)
fn central_sect2geo(elli: &Pelipsoide, phi1: f64, l1: f64, alpha: f64, dist: f64) -> (f64, f64) {
    let (sinl1, cosl1) = sincos(l1);
    let (sinphi1, cosphi1) = sincos(phi1);
    let beta1 = phi2beta(elli, phi1);
    let (sinbeta1, cosbeta1) = sincos(beta1);
    let (sinalpha, cosalpha) = sincos(alpha);

    let nx = (1.0 - elli.f) * (sinphi1 * sinl1 * cosalpha - cosl1 * sinalpha) * sinbeta1 + cosbeta1 * cosphi1 * sinl1 * cosalpha;
    let ny = (1.0 - elli.f) * (-sinphi1 * cosl1 * cosalpha - sinl1 * sinalpha) * sinbeta1 - cosbeta1 * cosphi1 * cosl1 * cosalpha;
    let nz = cosbeta1 * sinalpha;

    let psi0 = arg(nz, nx.hypot(ny));
    let l0 = if beta1 < 0.0 { arg(ny, -nx) } else { arg(-ny, nx) };

    let bdot = elli.a / (1.0 + elli.ep2 * psi0.sin().powi(2)).sqrt();
    let fdot = 1.0 - bdot / elli.a;
    let ellidot = Pelipsoide::new(fdot, elli.a);
    let platdot = Platn6::new(ellidot.n);
    let mudot = dist / ellidot.rmu;

    let sin_ll10 = (l1 - l0).sin();
    let cos_ll10 = (l1 - l0).cos();

    let u = ((1.0 - elli.e2) * sinbeta1.powi(2) + (cosbeta1 * sin_ll10).powi(2)).sqrt();
    let v = cosbeta1 * cos_ll10;
    let thetadot1 = arg(v, u);
    let mudot1 = lat2latn6(&platdot.geoc2rect, thetadot1);
    let mudot2 = mudot + mudot1;
    let thetadot2 = lat2latn6(&platdot.rect2geoc, mudot2);
    let thetadot = thetadot2 - thetadot1;

    let (sinthetadot, costhetadot) = sincos(thetadot);

    let (sinpsi1, cospsi1) = sincos(phi2psi(elli, phi1));
    let ux = cospsi1 * cosl1;
    let uy = cospsi1 * sinl1;
    let uz = sinpsi1;

    let vx = ny * uz - nz * uy;
    let vy = -nx * uz + nz * ux;
    let vz = nx * uy - ny * ux;
    let vnorma = (vx * vx + vy * vy + vz * vz).sqrt();
    let (wx, wy, wz) = (vx / vnorma, vy / vnorma, vz / vnorma);

    let rx = ux * costhetadot + wx * sinthetadot;
    let ry = uy * costhetadot + wy * sinthetadot;
    let rz = uz * costhetadot + wz * sinthetadot;

    let phi = arg((1.0 - elli.e2) * rx.hypot(ry), rz);
    let l = arg(rx, ry);

    (phi, l)
}

// ============================================================
// FUNCIONES INVERSAS (DISTANCIA Y ÁREA)
fn inv_align_dist_area(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_align(elli, beta1, l1, beta2, l2);
    let alpha1 = acimut_align(elli, &cs, beta1, l1);
    let salpha_limit = (PI / 12.0).sin().abs();
    let dist;
    let area;
    let phi0;
    let l0;
    let pathpoints: Vec<[f64; 5]>;
    if alpha1.sin().abs() < salpha_limit && (l2 - l1).abs() < PI_2 {
        (_, dist, area, phi0, l0, pathpoints) = inv_align_dist_area_beta(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    } else {
        (_, dist, area, phi0, l0, pathpoints) = inv_align_dist_area_lambda(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    }

    // Retorno final de la función
    (alpha1, dist, area, phi0, l0, pathpoints)
}

fn inv_align_dist_area_lambda(plat: &Platn6, elli: &Pelipsoide, alpha_ini: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;

    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_align(elli, beta1, l1, beta2, l2);
    let (beta0, l0) = calc_beta0_align(&cs, beta2, 0.5 * (l2 - l1));

    // Diferencia de longitud corregida
    let delta_l_total = tosemicirc(l2 - l1);
    let sgn_l = delta_l_total.signum();
    let delta_l_total_abs = delta_l_total.abs();

    let mut y = vec![beta1, 0.0, 0.0]; // [beta, s, S]
    let mut l = l1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_l_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();
    // Punto inicial
    //let alpha_ini = acimut_align(elli, &cs, beta11, l11);
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha_ini), 0.0, 0.0]);

    let derivs = |y: &[f64], l_val: f64| -> Vec<f64> {
        let beta = y[0];
        let c2beta = (2.0 * beta).cos();
        let (sbeta, cbeta) = sincos(beta);
        let (sl, cl) = sincos(l_val);

        let g_beta = c2beta * (cs[1] * cl + cs[2] * sl) + sbeta * (-cs[3] * cl - cs[4] * sl) + cs[5] * cbeta;
        let g_l = cbeta * (sbeta * (cs[2] * cl - cs[1] * sl) - cs[3] * sl + cs[4] * cl);
        let dbeta_dl = -g_l / (g_beta + sgn(g_beta) * EPSILON_MAQ);
        let ds_dl = elli.a * ((1.0 - elli.e2 * cbeta * cbeta) * dbeta_dl * dbeta_dl + cbeta * cbeta).sqrt();
        let tau = lat2latn6(&plat.para2auth, beta);
        let ds_dl_area = elli.r2 * elli.r2 * tau.sin();
        vec![dbeta_dl, ds_dl, ds_dl_area]
    };

    let mut iter_count = 0;
    while avanzado < delta_l_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_l_total_abs, h
            );
            break;
        }

        if avanzado + h > delta_l_total_abs {
            h = delta_l_total_abs - avanzado;
        }

        let h_signed = sgn_l * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, l);

        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado
            y = y_nuevo;
            l = tosemicirc(l + h_signed);
            avanzado += h;

            let alpha_val = acimut_align(elli, &cs, y[0], l);
            pathpoints.push([rad2deg(beta2phi(elli, y[0])), rad2deg(l), rad2deg(alpha_val), y[1].abs(), y[2]]);
        }
        // Si no se acepta, h se reduce abajo y se repite sin avanzar

        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    let alpha = acimut_align(elli, &cs, beta1, l1);
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);
    (alpha, dist, area, phi0, l0, pathpoints)
}

fn inv_align_dist_area_beta(plat: &Platn6, elli: &Pelipsoide, alpha_ini: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;
    // Transformaciones iniciales
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_align(elli, beta1, l1, beta2, l2);
    let (beta0, l0) = calc_beta0_align(&cs, beta2, 0.5 * (l2 - l1));

    // Diferencia de beta corregida y signo de integración
    let delta_beta_total = toquadcirc(beta2 - beta1);
    let sgn_beta = delta_beta_total.signum();

    // Caso extremo: misma latitud reducida (paralelo)
    let delta_beta_total_abs = delta_beta_total.abs();

    // Estado: y = [lambda, s, S] con beta como variable independiente
    let mut y = vec![l1, 0.0, 0.0];
    let mut beta = beta1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_beta_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();

    // Punto inicial
    //let alpha_ini = acimut_align(elli, &cs, beta11, l11);
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha_ini), 0.0, 0.0]);

    // Función de derivadas: dy/dbeta = f(beta, y)
    let derivs = |y: &[f64], beta_val: f64| -> Vec<f64> {
        let lambda = y[0];
        let (sbeta, cbeta) = sincos(beta_val);
        let (sl, cl) = sincos(lambda);

        // Cálculo de g_beta y g_lambda desde coeficientes de alineación
        let c2beta = (2.0 * beta_val).cos();
        let g_beta = c2beta * (cs[1] * cl + cs[2] * sl) + sbeta * (-cs[3] * cl - cs[4] * sl) + cs[5] * cbeta;
        let g_lambda = cbeta * (sbeta * (cs[2] * cl - cs[1] * sl) - cs[3] * sl + cs[4] * cl);

        // Derivada d(lambda)/d(beta) = -g_beta / g_lambda
        let dlambda_dbeta = -g_beta / (g_lambda + sgn(g_lambda) * EPSILON_MAQ);

        // Derivada de distancia: dL/dbeta = a * sqrt(1 - e²cos²β + cos²β·(dλ/dβ)²)
        let term = 1.0 - elli.e2 * cbeta * cbeta + cbeta * cbeta * dlambda_dbeta * dlambda_dbeta;
        let dL_dbeta = elli.a * term.sqrt();

        // Derivada de área: dS/dbeta = R_τ² · sin(τ(β)) · dλ/dβ
        let tau = lat2latn6(&plat.para2auth, beta_val);
        let dS_dbeta = elli.r2 * elli.r2 * tau.sin() * dlambda_dbeta;

        vec![dlambda_dbeta, dL_dbeta, dS_dbeta]
    };

    // Bucle de integración RK45 adaptativa
    let mut iter_count = 0;
    while avanzado < delta_beta_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align-beta) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_beta_total_abs, h
            );
            break;
        }

        // Ajuste del paso final para no sobrepasar el límite
        if avanzado + h > delta_beta_total_abs {
            h = delta_beta_total_abs - avanzado;
        }

        let h_signed = sgn_beta * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, beta);

        // Criterio de aceptación del paso (error relativo normalizado)
        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado: actualizar estado
            y = y_nuevo;
            beta = toquadcirc(beta + h_signed);
            avanzado += h;

            // Calcular acimut y guardar punto de trayectoria
            let alpha_val = acimut_align(elli, &cs, beta, y[0]);
            pathpoints.push([
                rad2deg(beta2phi(elli, beta)), // phi
                rad2deg(y[0]),                 // lambda
                rad2deg(alpha_val),            // alpha
                y[1].abs(),                    // distancia acumulada
                y[2],                          // área acumulada
            ]);
        }

        // Control adaptativo del tamaño de paso (orden 5 del método RK45)
        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    // Resultados finales
    let alpha = acimut_align(elli, &cs, beta1, l1);
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);

    (alpha, dist, area, phi0, l0, pathpoints)
}

fn inv_central_dist_area(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_central(elli, beta1, l1, beta2, l2);
    let alpha1 = acimut_central(elli, &cs, beta1, l1);
    let salpha_limit = (PI / 12.0).sin().abs();
    let dist;
    let area;
    let phi0;
    let l0;
    let pathpoints: Vec<[f64; 5]>;
    if alpha1.sin().abs() < salpha_limit && (l2 - l1).abs() < PI_2 {
        (_, dist, area, phi0, l0, pathpoints) = inv_central_dist_area_beta(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    } else {
        (_, dist, area, phi0, l0, pathpoints) = inv_central_dist_area_lambda(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    }

    // Retorno final de la función
    (alpha1, dist, area, phi0, l0, pathpoints)
}

fn inv_central_dist_area_lambda(plat: &Platn6, elli: &Pelipsoide, alpha1: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;

    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_central(elli, beta1, l1, beta2, l2);
    let beta11 = l2beta_central(&cs, beta1, l1);
    let (beta0, l0) = calc_beta0_central(&cs, beta2, 0.5 * (l2 - l1));

    let delta_l_total = tosemicirc(l2 - l1);
    let sgn_l = delta_l_total.signum();
    if sgn_l == 0.0 {
        let alpha = acimut_central(elli, &cs, beta1, l1);
        return (alpha, 0.0, 0.0, beta2phi(elli, beta0), l0, vec![]);
    }
    let delta_l_total_abs = delta_l_total.abs();

    let mut y = vec![beta11, 0.0, 0.0];
    let mut l = l1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_l_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();

    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha1), 0.0, 0.0]);

    let derivs = |y: &[f64], l_val: f64| -> Vec<f64> {
        let beta = y[0];
        let (sbeta, cbeta) = sincos(beta);
        let (sl, cl) = sincos(l_val);

        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let dbeta_dl = -g_l / g_beta;
        let ds_dl = elli.a * ((1.0 - elli.e2 * cbeta * cbeta) * dbeta_dl * dbeta_dl + cbeta * cbeta).sqrt();
        let tau = lat2latn6(&plat.para2auth, beta);
        let ds_dl_area = elli.r2 * elli.r2 * tau.sin();
        vec![dbeta_dl, ds_dl, ds_dl_area]
    };

    let mut iter_count = 0;

    while avanzado < delta_l_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (central) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_l_total_abs, h
            );
            break;
        }

        if avanzado + h > delta_l_total_abs {
            h = delta_l_total_abs - avanzado;
        }

        let h_signed = sgn_l * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, l);

        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            y = y_nuevo;
            l = tosemicirc(l + h_signed);
            avanzado += h;

            let alpha_val = acimut_central(elli, &cs, y[0], l);
            pathpoints.push([rad2deg(beta2phi(elli, y[0])), rad2deg(l), rad2deg(alpha_val), y[1].abs(), y[2]]);
        }

        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);
    (alpha1, dist, area, phi0, l0, pathpoints)
}

fn inv_central_dist_area_beta(plat: &Platn6, elli: &Pelipsoide, alpha1: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;
    // Transformaciones iniciales
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_central(elli, beta1, l1, beta2, l2);
    let (beta0, l0) = calc_beta0_central(&cs, beta2, 0.5 * (l2 - l1));

    // Diferencia de beta corregida y signo de integración
    let delta_beta_total = toquadcirc(beta2 - beta1);
    let sgn_beta = delta_beta_total.signum();

    // Caso extremo: misma latitud reducida (paralelo)
    let delta_beta_total_abs = delta_beta_total.abs();

    // Estado: y = [lambda, s, S] con beta como variable independiente
    let mut y = vec![l1, 0.0, 0.0];
    let mut beta = beta1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_beta_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();

    // Punto inicial
    //let alpha1 = acimut_align(elli, &cs, beta11, l11);
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha1), 0.0, 0.0]);

    // Función de derivadas: dy/dbeta = f(beta, y)
    let derivs = |y: &[f64], beta_val: f64| -> Vec<f64> {
        let lambda = y[0];
        let (sbeta, cbeta) = sincos(beta_val);
        let (sl, cl) = sincos(lambda);

        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_lambda = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;

        // Derivada d(lambda)/d(beta) = -g_beta / g_lambda
        let dlambda_dbeta = -g_beta / (g_lambda + sgn(g_lambda) * EPSILON_MAQ);

        // Derivada de distancia: dL/dbeta = a * sqrt(1 - e²cos²β + cos²β·(dλ/dβ)²)
        let term = 1.0 - elli.e2 * cbeta * cbeta + cbeta * cbeta * dlambda_dbeta * dlambda_dbeta;
        let dL_dbeta = elli.a * term.sqrt();

        // Derivada de área: dS/dbeta = R_τ² · sin(τ(β)) · dλ/dβ
        let tau = lat2latn6(&plat.para2auth, beta_val);
        let dS_dbeta = elli.r2 * elli.r2 * tau.sin() * dlambda_dbeta;

        vec![dlambda_dbeta, dL_dbeta, dS_dbeta]
    };

    // Bucle de integración RK45 adaptativa
    let mut iter_count = 0;
    while avanzado < delta_beta_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align-beta) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_beta_total_abs, h
            );
            break;
        }

        // Ajuste del paso final para no sobrepasar el límite
        if avanzado + h > delta_beta_total_abs {
            h = delta_beta_total_abs - avanzado;
        }

        let h_signed = sgn_beta * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, beta);

        // Criterio de aceptación del paso (error relativo normalizado)
        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado: actualizar estado
            y = y_nuevo;
            beta = toquadcirc(beta + h_signed);
            avanzado += h;

            // Calcular acimut y guardar punto de trayectoria
            let alpha_val = acimut_central(elli, &cs, beta, y[0]);
            pathpoints.push([
                rad2deg(beta2phi(elli, beta)), // phi
                rad2deg(y[0]),                 // lambda
                rad2deg(alpha_val),            // alpha
                y[1].abs(),                    // distancia acumulada
                y[2],                          // área acumulada
            ]);
        }

        // Control adaptativo del tamaño de paso (orden 5 del método RK45)
        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    // Resultados finales
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);

    (alpha1, dist, area, phi0, l0, pathpoints)
}

fn inv_normal_dist_area(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_normal(elli, beta1, l1, beta2, l2);
    let alpha1 = acimut_normal(elli, &cs, beta1, l1);
    let salpha_limit = (PI / 12.0).sin().abs();
    let dist;
    let area;
    let phi0;
    let l0;
    let pathpoints: Vec<[f64; 5]>;
    if alpha1.sin().abs() < salpha_limit && (l2 - l1).abs() < PI_2 {
        (_, dist, area, phi0, l0, pathpoints) = inv_normal_dist_area_beta(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    } else {
        (_, dist, area, phi0, l0, pathpoints) = inv_normal_dist_area_lambda(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    }

    // Retorno final de la función
    (alpha1, dist, area, phi0, l0, pathpoints)
}

fn inv_normal_dist_area_lambda(plat: &Platn6, elli: &Pelipsoide, alpha1: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;

    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_normal(elli, beta1, l1, beta2, l2);
    let (beta0, l0) = calc_beta0_normal(&cs, beta2, 0.5 * (l2 - l1));

    let delta_l_total = tosemicirc(l2 - l1);
    let sgn_l = delta_l_total.signum();

    let delta_l_total_abs = delta_l_total.abs();

    let mut y = vec![beta1, 0.0, 0.0];
    let mut l = l1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_l_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha1), 0.0, 0.0]);

    let derivs = |y: &[f64], l_val: f64| -> Vec<f64> {
        let beta = y[0];
        let (sbeta, cbeta) = sincos(beta);
        let (sl, cl) = sincos(l_val);

        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_l = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        let dbeta_dl = -g_l / g_beta;
        let ds_dl = elli.a * ((1.0 - elli.e2 * cbeta * cbeta) * dbeta_dl * dbeta_dl + cbeta * cbeta).sqrt();
        let tau = lat2latn6(&plat.para2auth, beta);
        let ds_dl_area = elli.r2 * elli.r2 * tau.sin();
        vec![dbeta_dl, ds_dl, ds_dl_area]
    };

    let mut iter_count = 0;

    while avanzado < delta_l_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (normal) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_l_total_abs, h
            );
            break;
        }

        if avanzado + h > delta_l_total_abs {
            h = delta_l_total_abs - avanzado;
        }

        let h_signed = sgn_l * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, l);

        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            y = y_nuevo;
            l = tosemicirc(l + h_signed);
            avanzado += h;

            let alpha_val = acimut_normal(elli, &cs, y[0], l);
            pathpoints.push([rad2deg(beta2phi(elli, y[0])), rad2deg(l), rad2deg(alpha_val), y[1].abs(), y[2]]);
        }

        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    let alpha = acimut_normal(elli, &cs, beta1, l1);
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);
    (alpha, dist, area, phi0, l0, pathpoints)
}

fn inv_normal_dist_area_beta(plat: &Platn6, elli: &Pelipsoide, alpha1: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;
    // Transformaciones iniciales
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let cs = cs_normal(elli, beta1, l1, beta2, l2);
    let (beta0, l0) = calc_beta0_normal(&cs, beta2, 0.5 * (l2 - l1));

    // Diferencia de beta corregida y signo de integración
    let delta_beta_total = toquadcirc(beta2 - beta1);
    let sgn_beta = delta_beta_total.signum();

    // Caso extremo: misma latitud reducida (paralelo)
    let delta_beta_total_abs = delta_beta_total.abs();

    // Estado: y = [lambda, s, S] con beta como variable independiente
    let mut y = vec![l1, 0.0, 0.0];
    let mut beta = beta1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_beta_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();

    // Punto inicial
    //let alpha1 = acimut_align(elli, &cs, beta11, l11);
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha1), 0.0, 0.0]);

    // Función de derivadas: dy/dbeta = f(beta, y)
    let derivs = |y: &[f64], beta_val: f64| -> Vec<f64> {
        let lambda = y[0];
        let (sbeta, cbeta) = sincos(beta_val);
        let (sl, cl) = sincos(lambda);

        let g_beta = -cs[1] * sbeta * cl - cs[2] * sbeta * sl + cs[3] * cbeta;
        let g_lambda = -cs[1] * cbeta * sl + cs[2] * cbeta * cl;
        // Derivada d(lambda)/d(beta) = -g_beta / g_lambda
        let dlambda_dbeta = -g_beta / (g_lambda + sgn(g_lambda) * EPSILON_MAQ);

        // Derivada de distancia: dL/dbeta = a * sqrt(1 - e²cos²β + cos²β·(dλ/dβ)²)
        let term = 1.0 - elli.e2 * cbeta * cbeta + cbeta * cbeta * dlambda_dbeta * dlambda_dbeta;
        let dL_dbeta = elli.a * term.sqrt();

        // Derivada de área: dS/dbeta = R_τ² · sin(τ(β)) · dλ/dβ
        let tau = lat2latn6(&plat.para2auth, beta_val);
        let dS_dbeta = elli.r2 * elli.r2 * tau.sin() * dlambda_dbeta;

        vec![dlambda_dbeta, dL_dbeta, dS_dbeta]
    };

    // Bucle de integración RK45 adaptativa
    let mut iter_count = 0;
    while avanzado < delta_beta_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align-beta) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_beta_total_abs, h
            );
            break;
        }

        // Ajuste del paso final para no sobrepasar el límite
        if avanzado + h > delta_beta_total_abs {
            h = delta_beta_total_abs - avanzado;
        }

        let h_signed = sgn_beta * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, beta);

        // Criterio de aceptación del paso (error relativo normalizado)
        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado: actualizar estado
            y = y_nuevo;
            beta = toquadcirc(beta + h_signed);
            avanzado += h;

            // Calcular acimut y guardar punto de trayectoria
            let alpha_val = acimut_normal(elli, &cs, beta, y[0]);
            pathpoints.push([
                rad2deg(beta2phi(elli, beta)), // phi
                rad2deg(y[0]),                 // lambda
                rad2deg(alpha_val),            // alpha
                y[1].abs(),                    // distancia acumulada
                y[2],                          // área acumulada
            ]);
        }

        // Control adaptativo del tamaño de paso (orden 5 del método RK45)
        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    // Resultados finales
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);

    (alpha1, dist, area, phi0, l0, pathpoints)
}

// Función principal para calcular la longitud del vértice en la línea geodésica, parte del código viene de QWEN.ai, debe ser revisado
fn calcular_lambda_vertice(elli: &Pelipsoide, phi1: f64, L1: f64, alpha1: f64, alpha0: f64) -> f64 {
    // 1. Obtener latitud paramétrica (beta) del punto 1
    let beta1 = phi2beta(elli, phi1);
    let (sbeta1, cbeta1) = beta1.sin_cos();

    // 2. Calcular sigma1 (arco esférico desde el ecuador hasta el punto 1)
    // tan(sigma1) = tan(beta1) / cos(alpha1)
    let sigma1 = sbeta1.atan2(alpha1.cos() * cbeta1);

    // 3. Determinar sigma en el vértice
    // Si el azimut va hacia el Norte (cos > 0), el vértice está a +PI/2
    // Si va hacia el Sur (cos < 0), está a -PI/2
    let sgn_cos_alpha1 = if alpha1.cos() >= 0.0 { 1.0 } else { -1.0 };
    let sigma_vertex = sgn_cos_alpha1 * PI_2;

    // Diferencia de arcos esféricos y sigma medio
    let sigma12_vertex = sigma_vertex - sigma1;
    let sigma_m_vertex = (sigma_vertex + sigma1) / 2.0;

    // 4. Calcular la longitud esférica (omega) para el punto 1 y el vértice
    // tan(omega) = sin(alpha0) * tan(sigma)
    let omega1 = (alpha0.sin() * sigma1.sin()).atan2(sigma1.cos());

    // En el vértice, tan(sigma) es infinito, por lo que omega = +/- PI/2
    let sgn_sin_alpha0 = if alpha0.sin() >= 0.0 { 1.0 } else { -1.0 };
    let omega_vertex = sgn_sin_alpha0 * PI_2;

    let delta_omega = omega_vertex - omega1;

    // 5. Calcular la corrección elipsoidal I para este tramo
    let n = elli.n;
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n2 * n2;
    let n5 = n4 * n;
    let n6 = n4 * n2;
    let c02 = sqr((alpha0).cos());
    let c04 = c02 * c02;
    let c06 = c04 * c02;
    let c08 = c04 * c04;
    let c010 = c08 * c02;

    // Coeficientes B (extraídos de tu L2omega12)
    let B0 = (-2.0 * n
        + (n2 + n3 + n4 + n5 + n6) * c02
        + (-(3.0 * n3) / 2.0 - (15.0 * n4) / 4.0 - (27.0 * n5) / 4.0 - (21.0 * n6) / 2.0) * c04
        + ((25.0 * n4) / 8.0 + (105.0 * n5) / 8.0 + 35.0 * n6) * c06
        + (-(245.0 * n5) / 32.0 - (735.0 * n6) / 16.0) * c08
        + ((1323.0 * n6) / 64.0) * c010)
        / (1.0 + n);

    let B1 = (n / 4.0 + n2 / 4.0 + n3 / 4.0 + n4 / 4.0 + n5 / 4.0) * c02
        + (-(3.0 * n2) / 8.0 - n3 - (15.0 * n4) / 8.0 - 3.0 * n5 + (5.0 * n6) / 8.0) * c04
        + ((51.0 * n3) / 64.0 + (229.0 * n4) / 64.0 + 10.0 * n5 - (95.0 * n6) / 16.0) * c06
        + (-(255.0 * n4) / 128.0 - (813.0 * n5) / 64.0 + (1161.0 * n6) / 64.0) * c08
        + ((701.0 * n5) / 128.0 - (11391.0 * n6) / 512.0) * c010;

    let B2 = (n2 / 16.0 + (5.0 * n3) / 32.0 + (9.0 * n4) / 32.0 + (7.0 * n5) / 16.0) * c04
        + (-(13.0 * n3) / 64.0 - (7.0 * n4) / 8.0 - (19.0 * n5) / 8.0 + (15.0 * n6) / 32.0) * c06
        + ((79.0 * n4) / 128.0 + (489.0 * n5) / 128.0 - (625.0 * n6) / 256.0) * c08
        + (-(487.0 * n5) / 256.0 + (1001.0 * n6) / 256.0) * c010;

    let B3 = ((5.0 * n3) / 192.0 + (7.0 * n4) / 64.0 + (7.0 * n5) / 24.0) * c06 + (-(17.0 * n4) / 128.0 - (155.0 * n5) / 192.0 + (41.0 * n6) / 192.0) * c08 + ((271.0 * n5) / 512.0 - (923.0 * n6) / 1536.0) * c010;

    let B4 = ((7.0 * n4) / 512.0 + (21.0 * n5) / 256.0) * c08 + (-(49.0 * n5) / 512.0 + (49.0 * n6) / 1024.0) * c010;

    let B5 = (21.0 * n5) / 2560.0 * c010;

    // Cálculo de S e I
    let S = B1 * (2.0 * sigma_m_vertex).cos() * (1.0 * sigma12_vertex).sin()
        + B2 * (4.0 * sigma_m_vertex).cos() * (2.0 * sigma12_vertex).sin()
        + B3 * (6.0 * sigma_m_vertex).cos() * (3.0 * sigma12_vertex).sin()
        + B4 * (8.0 * sigma_m_vertex).cos() * (4.0 * sigma12_vertex).sin()
        + B5 * (10.0 * sigma_m_vertex).cos() * (5.0 * sigma12_vertex).sin();

    let I_vertex = B0 * (sigma12_vertex + 2.0 * S);

    // 6. Convertir diferencia esférica a diferencia elipsoidal
    // Despejando L12 de tu fórmula: omega12 = L12 - alpha0.sin() * I
    let delta_L = delta_omega + alpha0.sin() * I_vertex;

    // 7. Longitud final del vértice
    let lambda0 = tosemicirc(L1 + delta_L);

    lambda0
}

fn inv_geodesic_dist_area_beta(plat: &Platn6, elli: &Pelipsoide, alpha1: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;

    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);

    // La nueva variable independiente es beta
    let delta_beta_total = beta2 - beta1;
    let sgn_beta = delta_beta_total.signum();
    let delta_beta_total_abs = delta_beta_total.abs();

    // Parámetros iniciales (seguimos usando l1, l2 para calcular el azimut inicial)
    let re0 = hypot(alpha1.cos(), alpha1.sin() * beta1.sin());
    let im0 = alpha1.sin() * beta1.cos();
    let alpha0 = arg(re0, im0);
    let lambda0 = calcular_lambda_vertice(&elli, phi1, l1, alpha1, alpha0);
    let beta0 = alpha0.sin().acos();

    // IMPORTANTE: El vector de estado 'y' ahora contiene [l, alpha, s, S].
    // Beta se ha convertido en la variable independiente (el "tiempo" del integrador).
    let mut y = vec![l1, alpha1, 0.0, 0.0];
    let mut beta = beta1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_beta_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha1), 0.0, 0.0]);

    let derivs = |y: &[f64], beta: f64| -> Vec<f64> {
        let alpha = y[1];

        let (sbeta, cbeta) = sincos(beta);
        let (salpha, calpha) = sincos(alpha);
        let D = (1.0 - elli.e2 * cbeta * cbeta).sqrt();

        // Derivadas originales respecto a l:
        // dbeta_dl = cbeta * calpha / salpha / D;
        // dalpha_dl = sbeta / D;
        // ds_dl = elli.a * cbeta / salpha;
        // dS_dl = elli.r2 * elli.r2 * tau.sin();

        // Nuevas derivadas respecto a beta (dx/dbeta = dx/dl / dbeta/dl):
        // ADVERTENCIA: En el vértice de la geodésica (alpha = pi/2), cos(alpha) = 0
        // y las derivadas respecto a beta tienden a infinito. Este método asume
        // que beta es monótona y que la geodésica NO cruza su vértice en este tramo.
        let tan_alpha = salpha / calpha;
        let tan_beta = sbeta / cbeta;

        let dl_dbeta = tan_alpha * D / cbeta;
        let dalpha_dbeta = tan_beta * tan_alpha;
        let ds_dbeta = elli.a * D / calpha;

        let tau = lat2latn6(&plat.para2auth, beta);
        let dS_dbeta = elli.r2 * elli.r2 * tau.sin() * dl_dbeta;

        vec![dl_dbeta, dalpha_dbeta, ds_dbeta, dS_dbeta]
    };

    let mut iter_count = 0;
    while avanzado < delta_beta_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_beta_total_abs, h
            );
            break;
        }

        if avanzado + h > delta_beta_total_abs {
            h = delta_beta_total_abs - avanzado;
        }

        let h_signed = sgn_beta * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, beta);

        let idist = MAXDIST.min((y_nuevo[2] - y[2]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado
            y = y_nuevo;
            beta += h_signed; // La latitud no necesita "tosemicirc" (no da la vuelta al elipsoide)
            avanzado += h;

            // NOTA: Aplicamos tosemicirc a y[0] (longitud) SOLO al guardar o imprimir.
            // Si lo aplicáramos dentro del vector 'y', crearíamos un salto de 2*PI
            // que rompería por completo el control de error adaptativo del integrador RK45.
            let l_wrapped = tosemicirc(y[0]);
            pathpoints.push([rad2deg(beta2phi(elli, beta)), rad2deg(l_wrapped), rad2deg(y[1]), y[2].abs(), y[3]]);
        }

        // Si no se acepta, h se reduce abajo y se repite sin avanzar

        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    let l_final = tosemicirc(y[0]);
    let dist = y[2];
    let area = y[3];
    println!("phi final = {:.12}", rad2deg(beta2phi(elli, beta)));
    println!("phi objetivo = {:.12}", rad2deg(phi2));
    println!("error phi = {:.6} m", (beta2phi(elli, beta) - phi2).abs() * elli.a);
    println!("L final = {:.12}", rad2deg(l_final));
    println!("L objetivo = {:.12}", rad2deg(l2));

    // Mejora: Usamos tosemicirc para calcular la diferencia angular más corta correctamente,
    // evitando errores de bulto si la ruta cruza el meridiano de +/- 180 grados.
    let err_l = tosemicirc(l_final - l2).abs();
    println!("error L = {:.6} m", err_l * elli.a);

    (alpha1, dist, area, beta0, lambda0, pathpoints)
}

fn inv_geodesic_dist_area_lambda(plat: &Platn6, elli: &Pelipsoide, alpha1: f64, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;

    let beta1 = phi2beta(elli, phi1);
    // Diferencia de longitud corregida
    let delta_l_total = tosemicirc(l2 - l1);
    let sgn_l = delta_l_total.signum();
    let re0 = hypot(alpha1.cos(), alpha1.sin() * beta1.sin());
    let im0 = alpha1.sin() * beta1.cos();
    let alpha0 = arg(re0, im0);
    let lambda0 = calcular_lambda_vertice(&elli, phi1, l1, alpha1, alpha0);
    let beta0 = alpha0.sin().acos();
    let delta_l_total_abs = delta_l_total.abs();

    let mut y = vec![beta1, alpha1, 0.0, 0.0]; // [beta, s, S]
    let mut l = l1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_l_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha1), 0.0, 0.0]);

    let derivs = |y: &[f64], _s: f64| -> Vec<f64> {
        let beta = y[0];
        let alpha = y[1];

        let (sbeta, cbeta) = sincos(beta);
        let (salpha, calpha) = sincos(alpha);
        let D = (1.0 - elli.e2 * sqr(cbeta)).sqrt();

        let dbeta_dl = cbeta * calpha / salpha / D;
        let dalpha_dl = sbeta / D;
        let ds_dl = elli.a * cbeta / salpha;
        let tau = lat2latn6(&plat.para2auth, beta);
        let dS_dl = elli.r2 * elli.r2 * tau.sin();

        vec![dbeta_dl, dalpha_dl, ds_dl, dS_dl]
    };
    let mut iter_count = 0;
    while avanzado < delta_l_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_l_total_abs, h
            );
            break;
        }

        if avanzado + h > delta_l_total_abs {
            h = delta_l_total_abs - avanzado;
        }

        let h_signed = sgn_l * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, l);

        let idist = MAXDIST.min((y_nuevo[2] - y[2]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado
            y = y_nuevo;
            l = tosemicirc(l + h_signed);
            avanzado += h;
            pathpoints.push([rad2deg(beta2phi(elli, y[0])), rad2deg(l), rad2deg(y[1]), y[2].abs(), y[3]]);
        }
        // Si no se acepta, h se reduce abajo y se repite sin avanzar

        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }

    let dist = y[2];
    let area = y[3];
    println!("phi final = {:.12}", rad2deg(beta2phi(elli, y[0])));
    println!("phi objetivo = {:.12}", rad2deg(phi2));
    println!("error phi = {:.6} m", (beta2phi(elli, y[0]) - phi2).abs() * elli.a);
    println!("L final = {:.12}", rad2deg(l));
    println!("L objetivo = {:.12}", rad2deg(l2));
    println!("error L = {:.6} m", (l - l2).abs() * elli.a);
    (alpha1, dist, area, beta0, lambda0, pathpoints)
}

fn inv_geodesic_dist_area(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let salpha_limit = (PI / 12.0).sin().abs();
    let phi0;
    let l0;
    let pathpoints: Vec<[f64; 5]>;
    let g = Geodesic::new(elli.a, elli.f);
    //let capabilities = caps::DISTANCE | caps::AZIMUTH | caps::REDUCEDLENGTH;
    let capabilities = caps::ALL;
    let (_, s12, mut alpha1, _azi2, _m12, _M12, _M21, mut area) = g._gen_inverse_azi(rad2deg(phi1), rad2deg(l1), rad2deg(phi2), rad2deg(l2), capabilities);
    println!("area geodesica: {}", area);
    alpha1 = deg2rad(alpha1);
    if alpha1.sin().abs() < salpha_limit && (l2 - l1).abs() < PI_2 {
        //(alpha1, dist, area, phi0, l0, pathpoints)=
        (_, _, area, phi0, l0, pathpoints) = inv_geodesic_dist_area_beta(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    } else {
        (_, _, area, phi0, l0, pathpoints) = inv_geodesic_dist_area_lambda(plat, elli, alpha1, phi1, l1, phi2, l2, max_step_ang);
    }
    (alpha1, s12, area, phi0, l0, pathpoints)
}

fn inv_loxo_dist_area(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let chi1 = lat2latn6(&plat.para2conf, beta1);
    let chi2 = lat2latn6(&plat.para2conf, beta2);
    let q1 = lam(chi1);
    let q2 = lam(chi2);
    let alpha = arg(q2 - q1, tosemicirc(l2 - l1));
    let salpha_limit = (PI / 12.0).sin().abs();
    let dist;
    let area;
    let phi0;
    let l0;
    let pathpoints: Vec<[f64; 5]>;
    if alpha.sin().abs() < salpha_limit && (l2 - l1).abs() < PI_2 {
        (_, dist, area, phi0, l0, pathpoints) = inv_loxo_dist_area_beta(plat, elli, alpha, phi1, l1, phi2, l2, max_step_ang);
    } else {
        (_, dist, area, phi0, l0, pathpoints) = inv_loxo_dist_area_lambda(plat, elli, alpha, phi1, l1, phi2, l2, max_step_ang);
    }
    (alpha, dist, area, phi0, l0, pathpoints)
}

fn inv_loxo_dist_area_lambda(plat: &Platn6, elli: &Pelipsoide, alpha: f64, phi1: f64, l1: f64, _phi2: f64, l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;

    let (salpha, calpha) = sincos(alpha);
    let beta1 = phi2beta(elli, phi1);
    let (phi0, l0) = (PI_2, 0.0);

    // Diferencia de longitud corregida
    let delta_l_total = tosemicirc(l2 - l1);
    let sgn_l = delta_l_total.signum();
    let delta_l_total_abs = delta_l_total.abs();

    let mut y = vec![beta1, 0.0, 0.0]; // [beta, s, S]
    let mut l = l1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_l_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha), 0.0, 0.0]);

    let derivs = |y: &[f64], _l_val: f64| -> Vec<f64> {
        let beta = y[0];
        let (_sbeta, cbeta) = sincos(beta);
        let D = (1.0 - elli.e2 * cbeta * cbeta).sqrt();

        let dbeta_dl = div(cbeta * calpha, D * salpha);
        let ds_dl = elli.a * (D * D * dbeta_dl.powi(2) + cbeta * cbeta).sqrt();
        let tau = lat2latn6(&plat.para2auth, beta);
        let ds_dl_area = elli.r2 * elli.r2 * tau.sin();
        vec![dbeta_dl, ds_dl, ds_dl_area]
    };

    let mut iter_count = 0;
    while avanzado < delta_l_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_l_total_abs, h
            );
            break;
        }

        if avanzado + h > delta_l_total_abs {
            h = delta_l_total_abs - avanzado;
        }

        let h_signed = sgn_l * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, l);

        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);

        if err_norm < tol || h <= h_min {
            // Paso aceptado
            y = y_nuevo;
            l = tosemicirc(l + h_signed);
            avanzado += h;
            pathpoints.push([rad2deg(beta2phi(elli, y[0])), rad2deg(l), rad2deg(alpha), y[1].abs(), y[2]]);
        }
        // Si no se acepta, h se reduce abajo y se repite sin avanzar

        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }
    let dist = y[1].abs();
    let area = y[2];
    (alpha, dist, area, phi0, l0, pathpoints)
}

fn inv_loxo_dist_area_beta(plat: &Platn6, elli: &Pelipsoide, alpha: f64, phi1: f64, l1: f64, phi2: f64, _l2: f64, max_step_ang: f64) -> (f64, f64, f64, f64, f64, Vec<[f64; 5]>) {
    let tol = EPSILON_ANG;
    let h_max = max_step_ang;
    let h_min = EPSILON_ANG;
    // Transformaciones iniciales
    let beta1 = phi2beta(elli, phi1);
    let beta2 = phi2beta(elli, phi2);
    let (beta0, l0) = (PI_2, 0.0);

    // Diferencia de beta corregida y signo de integración
    let delta_beta_total = toquadcirc(beta2 - beta1);
    let sgn_beta = delta_beta_total.signum();

    // Caso extremo: misma latitud reducida (paralelo)
    let delta_beta_total_abs = delta_beta_total.abs();

    // Estado: y = [lambda, s, S] con beta como variable independiente
    let mut y = vec![l1, 0.0, 0.0];
    let mut beta = beta1;
    let mut avanzado = 0.0;
    let mut h = h_max.min(delta_beta_total_abs);

    let mut pathpoints: Vec<[f64; 5]> = Vec::new();
    pathpoints.push([rad2deg(beta2phi(elli, beta1)), rad2deg(l1), rad2deg(alpha), 0.0, 0.0]);
    // Función de derivadas: dy/dbeta = f(beta, y)
    let (salpha, calpha) = sincos(alpha);
    let derivs = |_y: &[f64], beta_val: f64| -> Vec<f64> {
        let (_sbeta, cbeta) = sincos(beta_val);
        let D = (1.0 - elli.e2 * cbeta * cbeta).sqrt();
        let dlambda_dbeta = div(D * salpha, cbeta * calpha);
        let ds_dbeta = elli.a * (D * D + cbeta * cbeta * dlambda_dbeta.powi(2)).sqrt();
        // Derivada de área: dS/dbeta = R_τ² · sin(τ(β)) · dλ/dβ
        let tau = lat2latn6(&plat.para2auth, beta_val);
        let dS_dbeta = elli.r2 * elli.r2 * tau.sin() * dlambda_dbeta;
        vec![dlambda_dbeta, ds_dbeta, dS_dbeta]
    };

    // Bucle de integración RK45 adaptativa
    let mut iter_count = 0;
    while avanzado < delta_beta_total_abs - EPSILON_ANG {
        iter_count += 1;
        if iter_count > MAX_ITER {
            eprintln!(
                "Error: Integración no converge (align-beta) tras {} iteraciones.\n  Avance: {:.6} de {:.6}, h actual: {:.2e}",
                MAX_ITER, avanzado, delta_beta_total_abs, h
            );
            break;
        }
        // Ajuste del paso final para no sobrepasar el límite
        if avanzado + h > delta_beta_total_abs {
            h = delta_beta_total_abs - avanzado;
        }

        let h_signed = sgn_beta * h;
        let (y_nuevo, error1, _error0) = rk45_step(&mut y, derivs, h_signed, beta);
        // Criterio de aceptación del paso (error relativo normalizado)
        let idist = MAXDIST.min((y_nuevo[1] - y[1]).abs());
        let err_norm = error1 / idist.max(EPSILON_ANG);
        if err_norm < tol || h <= h_min {
            // Paso aceptado: actualizar estado
            y = y_nuevo;
            beta = toquadcirc(beta + h_signed);
            avanzado += h;
            // Calcular acimut y guardar punto de trayectoria
            pathpoints.push([
                rad2deg(beta2phi(elli, beta)), // phi
                rad2deg(y[0]),                 // lambda
                rad2deg(alpha),                // alpha
                y[1].abs(),                    // distancia acumulada
                y[2],                          // área acumulada
            ]);
        }

        // Control adaptativo del tamaño de paso (orden 5 del método RK45)
        let factor = if err_norm == 0.0 { 5.0 } else { (tol / (err_norm + EPSILON_ANG)).powf(1.0 / 5.0) };
        h = h_max.min(h_min.max(h * factor));
    }
    let dist = y[1].abs();
    let area = y[2];
    let phi0 = beta2phi(elli, beta0);

    (alpha, dist, area, phi0, l0, pathpoints)
}

fn direct_curva_geodesic(_dummy_plat: &Platn6, elli: &Pelipsoide, phi1: f64, lambda1: f64, alpha1: f64, s12: f64, _dummy_max_step_ang: f64) -> (f64, f64) {
    let g = Geodesic::new(elli.a, elli.f);
    let (phi2, lambda2, _alpha1) = g.direct(rad2deg(phi1), rad2deg(lambda1), rad2deg(alpha1), s12);
    //let (_alpha0, _alpha1, phi2, lambda2) = karney_direct_curva_geodesic(elli, phi1, lambda1, alpha1, s12);
    (deg2rad(phi2), deg2rad(lambda2))
}

fn alpha_s_loxo(plat: &Platn6, elli: &Pelipsoide, beta1: f64, lambda1: f64, beta2: f64, lambda2: f64) -> (f64, f64) {
    let chi1 = lat2latn6(&plat.para2conf, beta1);
    let chi2 = lat2latn6(&plat.para2conf, beta2);
    let q1 = lam(chi1);
    let q2 = lam(chi2);
    let L12 = tosemicirc(lambda2 - lambda1);
    let alpha = arg(q2 - q1, L12);
    let (salpha, calpha) = sincos(alpha);
    let salpha_limit = (PI / 12.0).sin().abs();
    let dist: f64;
    if salpha.abs() < salpha_limit {
        let mu1 = lat2latn6(&plat.para2rect, beta1);
        let mu2 = lat2latn6(&plat.para2rect, beta2);
        let s1 = elli.rmu * mu1;
        let s2 = elli.rmu * mu2;
        dist = ((s2 - s1) / calpha).abs();
    } else {
        let cbeta1 = beta1.cos();
        let cbeta2 = beta2.cos();
        let s1 = elli.a * cbeta1 * L12;
        let s2 = elli.a * cbeta2 * L12;
        dist = ((s2 - s1) / salpha).abs();
    }
    (alpha, dist)
}

fn direct_curva_meridian(plat: &Platn6, elli: &Pelipsoide, phi1: f64, lambda1: f64, alpha: f64, s12: f64) -> (f64, f64) {
    let mu1 = lat2latn6(&plat.geod2rect, phi1);
    let mut lambda2 = lambda1;
    let sgn_s12 = sgn(alpha.sin());
    let s1 = elli.rmu * mu1;
    let Q = elli.rmu * PI_2;
    let s2 = s12 * sgn_s12 + s1;
    let mu2 = s2 / elli.rmu;
    let phi2 = lat2latn6(&plat.rect2geod, mu2);
    if s2.abs() > Q {
        lambda2 = tosemicirc(lambda1 + PI);
    }
    (phi2, lambda2)
}

fn dist_curva_meridian(plat: &Platn6, elli: &Pelipsoide, phi1: f64, lambda1: f64, phi2: f64, lambda2: f64) -> f64 {
    let mu1 = lat2latn6(&plat.geod2rect, phi1);
    let mu2 = lat2latn6(&plat.geod2rect, phi2);
    let s1 = mu1 * elli.rmu;
    let s2 = mu2 * elli.rmu;
    let s12;
    let Q = PI_2 * elli.rmu;
    let L12 = tosemicirc(lambda2 - lambda1);
    if L12.abs() > (PI - EPSILON_ANG) {
        s12 = Q * sgn(mu2) - mu2 - Q * sgn(mu1) - mu1;
    } else {
        s12 = s2 - s1;
    }
    s12
}

fn direct_curva_loxo(plat: &Platn6, elli: &Pelipsoide, phi1: f64, lambda1: f64, alpha: f64, s12: f64, _dummy_max_step_ang: f64) -> (f64, f64) {
    let beta1 = phi2beta(elli, phi1);
    let lambda2;
    let phi2;
    let (salpha, calpha) = sincos(alpha);

    if salpha.abs() < 1.0 {
        if salpha.abs() <= EPSILON_ANG {
            (phi2, lambda2) = direct_curva_meridian(plat, elli, phi1, lambda1, alpha, s12);
        } else {
            let mu1 = lat2latn6(&plat.para2rect, beta1);
            let s1 = elli.rmu * mu1;
            let mu2 = (s12 * calpha + s1) / elli.rmu;
            phi2 = lat2latn6(&plat.rect2geod, mu2);
            let chi1 = lat2latn6(&plat.para2conf, beta1);
            let chi2 = lat2latn6(&plat.geod2conf, phi2);
            let q1 = lam(chi1);
            let q2 = lam(chi2);
            lambda2 = tosemicirc(salpha / calpha * (q2 - q1) + lambda1);
        }
    } else {
        //la loxo es un paralelo
        let cbeta1 = beta1.cos();
        let r = elli.a * cbeta1;
        let L12 = s12 / r * sgn(alpha);
        lambda2 = tosemicirc(L12 + lambda1);
        phi2 = beta2phi(elli, beta1);
    }
    (phi2, lambda2)
}

fn direct_curva_align(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, alpha_target: f64, dist_target: f64, max_step_ang: f64) -> (f64, f64) {
    let mut phi2;
    let mut l2;
    if alpha_target.sin().abs() <= EPSILON_ANG {
        (phi2, l2) = direct_curva_meridian(plat, elli, phi1, l1, alpha_target, dist_target);
    } else {
        (phi2, l2) = central_sect2geo(elli, phi1, l1, alpha_target, dist_target);
        let mut beta2 = phi2beta(elli, phi2);
        for _ in 0..50 {
            // Calcular inverso con perturbaciones
            let (alpha, dist, _, _, _, _) = inv_align_dist_area(plat, elli, phi1, l1, phi2, l2, max_step_ang);
            let (alpha_b, dist_b, _, _, _, _) = inv_align_dist_area(plat, elli, phi1, l1, beta2phi(elli, beta2 + MAXSTEP_H), l2, max_step_ang);
            let (alpha_l, dist_l, _, _, _, _) = inv_align_dist_area(plat, elli, phi1, l1, phi2, l2 + MAXSTEP_H, max_step_ang);

            let a = (alpha_b - alpha) / MAXSTEP_H;
            let b = (alpha_l - alpha) / MAXSTEP_H;
            let c = (dist_b - dist) / MAXSTEP_H;
            let d = (dist_l - dist) / MAXSTEP_H;

            let f1 = alpha - alpha_target;
            let f2 = dist - dist_target;
            let det = a * d - b * c;
            let sgndet = sgn(det);

            let mut delta_beta = (d * f1 - b * f2) / (det + sgndet * EPSILON_MAQ);
            let mut delta_l = (-c * f1 + a * f2) / (det + sgndet * EPSILON_MAQ);

            if delta_l.abs() > MAXSTEP_ANG {
                delta_l = MAXSTEP_ANG.copysign(delta_l);
            }
            if delta_beta.abs() > MAXSTEP_ANG {
                delta_beta = MAXSTEP_ANG.copysign(delta_beta);
            }

            beta2 -= delta_beta;
            l2 -= delta_l;
            phi2 = beta2phi(elli, beta2);

            if delta_beta.abs().max(delta_l.abs()) < EPSILON_ANG {
                break;
            }
        }
    }
    (toquadcirc(phi2), tosemicirc(l2))
}

fn direct_curva_central(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, alpha_target: f64, dist_target: f64, max_step_ang: f64) -> (f64, f64) {
    let mut phi2;
    let mut l2;
    if alpha_target.sin().abs() <= EPSILON_ANG {
        (phi2, l2) = direct_curva_meridian(plat, elli, phi1, l1, alpha_target, dist_target);
    } else {
        (phi2, l2) = central_sect2geo(elli, phi1, l1, alpha_target, dist_target);
        //let beta1 = phi2beta(elli, phi1);
        let mut beta2 = phi2beta(elli, phi2);

        for _ in 0..50 {
            let (alpha, dist, _, _, _, _) = inv_central_dist_area(plat, elli, phi1, l1, phi2, l2, max_step_ang);
            let (alpha_b, dist_b, _, _, _, _) = inv_central_dist_area(plat, elli, phi1, l1, beta2phi(elli, beta2 + MAXSTEP_H), l2, max_step_ang);
            let (alpha_l, dist_l, _, _, _, _) = inv_central_dist_area(plat, elli, phi1, l1, phi2, l2 + MAXSTEP_H, max_step_ang);

            let a = (alpha_b - alpha) / MAXSTEP_H;
            let b = (alpha_l - alpha) / MAXSTEP_H;
            let c = (dist_b - dist) / MAXSTEP_H;
            let d = (dist_l - dist) / MAXSTEP_H;

            let f1 = alpha - alpha_target;
            let f2 = dist - dist_target;
            let det = a * d - b * c;
            let sgndet = sgn(det);

            let mut delta_beta = (d * f1 - b * f2) / (det + sgndet * EPSILON_MAQ);
            let mut delta_l = (-c * f1 + a * f2) / (det + sgndet * EPSILON_MAQ);

            if delta_l.abs() > MAXSTEP_ANG {
                delta_l = MAXSTEP_ANG.copysign(delta_l);
            }
            if delta_beta.abs() > MAXSTEP_ANG {
                delta_beta = MAXSTEP_ANG.copysign(delta_beta);
            }

            beta2 -= delta_beta;
            phi2 = beta2phi(elli, beta2);
            l2 -= delta_l;

            if delta_beta.abs().max(delta_l.abs()) < EPSILON_ANG {
                break;
            }
        }
    }
    (toquadcirc(phi2), tosemicirc(l2))
}

fn direct_curva_normal(plat: &Platn6, elli: &Pelipsoide, phi1: f64, l1: f64, alpha_target: f64, dist_target: f64, max_step_ang: f64) -> (f64, f64) {
    let mut phi2;
    let mut l2;
    if alpha_target.sin().abs() <= EPSILON_ANG {
        (phi2, l2) = direct_curva_meridian(plat, elli, phi1, l1, alpha_target, dist_target);
    } else {
        (phi2, l2) = central_sect2geo(elli, phi1, l1, alpha_target, dist_target);
        let mut beta2 = phi2beta(elli, phi2);

        for _ in 0..50 {
            let (alpha, dist, _, _, _, _) = inv_normal_dist_area(plat, elli, phi1, l1, phi2, l2, max_step_ang);
            let (alpha_b, dist_b, _, _, _, _) = inv_normal_dist_area(plat, elli, phi1, l1, beta2phi(elli, beta2 + MAXSTEP_H), l2, max_step_ang);
            let (alpha_l, dist_l, _, _, _, _) = inv_normal_dist_area(plat, elli, phi1, l1, phi2, l2 + MAXSTEP_H, max_step_ang);

            let a = (alpha_b - alpha) / MAXSTEP_H;
            let b = (alpha_l - alpha) / MAXSTEP_H;
            let c = (dist_b - dist) / MAXSTEP_H;
            let d = (dist_l - dist) / MAXSTEP_H;

            let f1 = alpha - alpha_target;
            let f2 = dist - dist_target;
            let det = a * d - b * c;
            let sgndet = sgn(det);

            let mut delta_beta = (d * f1 - b * f2) / (det + sgndet * EPSILON_MAQ);
            let mut delta_l = (-c * f1 + a * f2) / (det + sgndet * EPSILON_MAQ);

            if delta_l.abs() > MAXSTEP_ANG {
                delta_l = MAXSTEP_ANG.copysign(delta_l);
            }
            if delta_beta.abs() > MAXSTEP_ANG {
                delta_beta = MAXSTEP_ANG.copysign(delta_beta);
            }

            beta2 -= delta_beta;
            l2 -= delta_l;
            phi2 = beta2phi(elli, beta2);

            if delta_beta.abs().max(delta_l.abs()) < EPSILON_ANG {
                break;
            }
        }
    }
    (toquadcirc(phi2), tosemicirc(l2))
}

// ============================================================
// GUARDAR KMZ (SIMPLIFICADO)

fn guardar_kmz(pathpoints: &[[f64; 5]], filename: &str) -> io::Result<()> {
    let mut kml_content = String::from(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
<name>Curva Geodesica</name>
<Style id="redLine">
<LineStyle>
<color>ff0000ff</color>
<width>3</width>
</LineStyle>
</Style>
<Placemark>
<name>Curva</name>
<styleUrl>#redLine</styleUrl>
<LineString>
<coordinates>"#,
    );

    for point in pathpoints {
        kml_content.push_str(&format!("{},{},0.0\n", point[1], point[0]));
    }

    kml_content.push_str(
        r#"</coordinates>
</LineString>
</Placemark>
</Document>
</kml>"#,
    );

    let mut file = File::create(filename)?;
    file.write_all(kml_content.as_bytes())?;
    Ok(())
}

fn guardar_shp(pathpoints: &[[f64; 5]], filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    // ---------------------------------------------------------
    // Geometría
    // ---------------------------------------------------------

    let points: Vec<Point> = pathpoints
        .iter()
        .map(|p| Point::new(p[1], p[0])) // lon, lat
        .collect();

    let polyline = Polyline::new(points);

    // ---------------------------------------------------------
    // Tabla DBF
    // ---------------------------------------------------------

    let table_builder = TableWriterBuilder::new()
        .add_character_field(FieldName::try_from("TIPO")?, 20)
        .add_numeric_field(FieldName::try_from("NPTS")?, 10, 0);

    let mut writer = Writer::from_path(filename, table_builder)?;

    let mut record = Record::default();

    record.insert("TIPO".to_string(), FieldValue::Character(Some("CURVA".to_string())));

    record.insert("NPTS".to_string(), FieldValue::Numeric(Some(pathpoints.len() as f64)));

    writer.write_shape_and_record(&polyline, &record)?;

    // Asegurar escritura completa
    drop(writer);

    // ---------------------------------------------------------
    // Archivo PRJ (WGS84)
    // ---------------------------------------------------------
    /*
        let prj_filename = filename.replace(".shp", ".prj");

        let wgs84_wkt = r#"GEOGCS["WGS 84",
DATUM["WGS_1984",
SPHEROID["WGS 84",6378137,298.257223563]],
PRIMEM["Greenwich",0],
UNIT["degree",0.0174532925199433]]"#;

        let mut prj = File::create(prj_filename)?;
        prj.write_all(wgs84_wkt.as_bytes())?;
    */
    Ok(())
}

fn guardar_csv(pathpoints: &[[f64; 5]], filename: &str) -> io::Result<()> {
    let mut wtr = csv::Writer::from_path(filename)?;
    wtr.write_record(&["Latitud", "Longitud", "Acimut", "Distancia", "Area"])?;
    for point in pathpoints {
        wtr.write_record(&[
            format!("{:.10}", point[0]),
            format!("{:.10}", point[1]),
            format!("{:.10}", point[2]),
            format!("{:.4}", point[3]),
            format!("{:.4}", point[4]),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

fn guardar_shp_points(pathpoints: &[[f64; 5]], filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    // curva.shp -> curva_points.shp
    let shp_filename = if filename.ends_with(".shp") {
        filename.replace(".shp", "_points.shp")
    } else {
        format!("{filename}_points.shp")
    };

    let table_builder = TableWriterBuilder::new()
        .add_numeric_field(FieldName::try_from("LATITUD")?, 20, 10)
        .add_numeric_field(FieldName::try_from("LONGITUD")?, 20, 10)
        .add_numeric_field(FieldName::try_from("ACIMUT")?, 20, 10)
        .add_numeric_field(FieldName::try_from("DISTANCIA")?, 20, 4)
        .add_numeric_field(FieldName::try_from("AREA")?, 20, 4);

    let mut writer = Writer::from_path(&shp_filename, table_builder)?;

    for p in pathpoints {
        let point = Point::new(p[1], p[0]); // lon, lat
        let mut record = Record::default();
        record.insert("LATITUD".to_string(), FieldValue::Numeric(Some(p[0])));
        record.insert("LONGITUD".to_string(), FieldValue::Numeric(Some(p[1])));
        record.insert("ACIMUT".to_string(), FieldValue::Numeric(Some(p[2])));
        record.insert("DISTANCIA".to_string(), FieldValue::Numeric(Some(p[3])));
        record.insert("AREA".to_string(), FieldValue::Numeric(Some(p[4])));
        writer.write_shape_and_record(&point, &record)?;
    }

    drop(writer);
    /*
        // Crear .prj WGS84
        let prj_filename = Path::new(&shp_filename).with_extension("prj");

        let wkt = r#"GEOGCS["WGS 84",
DATUM["WGS_1984",
SPHEROID["WGS 84",6378137,298.257223563]],
PRIMEM["Greenwich",0],
UNIT["degree",0.0174532925199433]]"#;

        let mut prj = File::create(prj_filename)?;
        prj.write_all(wkt.as_bytes())?;
    */
    Ok(())
}

// ============================================================
// lectura de polígonos, siempre el orden es latitud, longitud en grados decimales

/// Lee un polígono desde un archivo CSV/TXT, detectando automáticamente el formato.
/// Soporta:
/// 1. `Lat, Lon` o `Lat, Lon, Alt`
/// 2. `ID, Lat, Lon, Alt` (donde ID es un entero)
pub fn leer_poligono(filepath: &str) -> io::Result<Vec<(f64, f64)>> {
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);
    let mut vertices = Vec::new();

    for (_line_num, line) in reader.lines().enumerate() {
        let line = line?;

        // 1. Limpiar espacios y posible BOM (Byte Order Mark) de archivos UTF-8
        let line = line.trim().trim_start_matches('\u{FEFF}');

        // Ignorar líneas vacías o comentarios
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // 2. Separar por comas, tabulaciones o punto y coma
        let parts: Vec<&str> = line.split(|c: char| c == ',' || c == '\t' || c == ';').map(|s| s.trim()).filter(|s| !s.is_empty()).collect();

        // 3. Intentar extraer las coordenadas usando la heurística
        if let Some((lat, lon)) = extraer_coordenadas(&parts) {
            vertices.push((lat, lon));
        }
        // Si no coincide con ningún formato, simplemente ignora la línea (ej. encabezados)
    }

    if vertices.len() < 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Se requieren al menos 3 vértices válidos para un polígono, pero se encontraron {}", vertices.len()),
        ));
    }

    Ok(vertices)
}

/// Función auxiliar que aplica la heurística para extraer (Lat, Lon)
fn extraer_coordenadas(parts: &[&str]) -> Option<(f64, f64)> {
    // Closures auxiliares para validación
    let es_probable_id = |s: &str| s.parse::<i64>().is_ok() && !s.contains('.');
    let es_lat_valida = |v: f64| (-90.0..=90.0).contains(&v);
    let es_lon_valida = |v: f64| (-180.0..=180.0).contains(&v);

    // CASO 1: Formato con ID entero (ID, Lat, Lon, [Alt])
    // Verificamos que la primera columna sea un entero para no confundirla con una Latitud
    if parts.len() >= 3 && es_probable_id(parts[0]) {
        if let (Ok(lat), Ok(lon)) = (parts[1].parse::<f64>(), parts[2].parse::<f64>()) {
            if es_lat_valida(lat) && es_lon_valida(lon) {
                return Some((lat, lon));
            }
        }
    }

    // CASO 2: Formato sin ID (Lat, Lon, [Alt])
    if parts.len() >= 2 {
        if let (Ok(lat), Ok(lon)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
            if es_lat_valida(lat) && es_lon_valida(lon) {
                return Some((lat, lon));
            }
        }
    }

    // Si no encaja en ningún caso (ej. fila de encabezados como "id,lat,lon,alt"), retorna None
    None
}

// ============================================================
// MAIN

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        print_help();
        process::exit(1);
    }

    let mut i = 1;
    let mut tipo = "central".to_string();
    let mut output_base = String::new();
    let mut max_step = MAXSTEPRK45_ANG_INITIALVALUE;
    let mut p1: Option<(f64, f64)> = None;
    let mut p2: Option<(f64, f64)> = None;
    let mut azimut: Option<f64> = None;
    let mut distancia: Option<f64> = None;
    let mut modo = String::new();
    let mut poly_file: Option<String> = None;
    let mut semi_major = GRS80_A;
    let mut inv_f = 298.2572221008827;

    while i < args.len() {
        let arg = args[i].to_lowercase();
        match arg.as_str() {
            "-i" | "--inverso" => {
                modo = "inverso".to_string();
                i += 1;
            }
            "-d" | "--directo" => {
                modo = "directo".to_string();
                i += 1;
            }
            "-poly" | "--poly-sup" => {
                modo = "poly".to_string();
                i += 1;
                if i < args.len() {
                    poly_file = Some(args[i].clone());
                    i += 1;
                } else {
                    eprintln!("Error: -poly requiere la ruta de un archivo");
                    process::exit(1);
                }
            }
            "-p1" => {
                i += 1;
                if i + 1 < args.len() {
                    if let (Ok(lat), Ok(lon)) = (args[i].parse::<f64>(), args[i + 1].parse::<f64>()) {
                        p1 = Some((lat, lon));
                        i += 2;
                    } else {
                        eprintln!("Error: -P1 necesita dos números (latitud, longitud)");
                        process::exit(1);
                    }
                } else {
                    eprintln!("Error: -P1 necesita dos números");
                    process::exit(1);
                }
            }
            "-p2" => {
                i += 1;
                if i + 1 < args.len() {
                    if let (Ok(lat), Ok(lon)) = (args[i].parse::<f64>(), args[i + 1].parse::<f64>()) {
                        p2 = Some((lat, lon));
                        i += 2;
                    } else {
                        eprintln!("Error: -P2 necesita dos números");
                        process::exit(1);
                    }
                } else {
                    eprintln!("Error: -P2 necesita dos números");
                    process::exit(1);
                }
            }
            "-az" | "--azimuth" => {
                i += 1;
                if i < args.len() {
                    if let Ok(az) = args[i].parse::<f64>() {
                        azimut = Some(az);
                        i += 1;
                    } else {
                        eprintln!("Error: -az necesita un número (acimut en grados)");
                        process::exit(1);
                    }
                } else {
                    eprintln!("Error: -az necesita un número");
                    process::exit(1);
                }
            }
            "-s" | "--distance" => {
                i += 1;
                if i < args.len() {
                    if let Ok(dist) = args[i].parse::<f64>() {
                        distancia = Some(dist);
                        i += 1;
                    } else {
                        eprintln!("Error: -s necesita un número (distancia en metros)");
                        process::exit(1);
                    }
                } else {
                    eprintln!("Error: -s necesita un número");
                    process::exit(1);
                }
            }
            "-t" | "--tipo" => {
                i += 1;
                if i < args.len() {
                    tipo = args[i].to_lowercase();
                    i += 1;
                    if !["align", "central", "normal", "geodesic", "loxo"].contains(&tipo.as_str()) {
                        eprintln!("Error: tipo debe ser align, central, geodesic o normal");
                        process::exit(1);
                    }
                }
            }
            "-e" | "--ellipsoid" => {
                i += 1;

                if i + 1 < args.len() {
                    if let (Ok(a), Ok(invf)) = (args[i].parse::<f64>(), args[i + 1].parse::<f64>()) {
                        semi_major = a;
                        inv_f = invf;
                        i += 2;
                    } else {
                        eprintln!("Error: -e necesita dos números: a inv_f");
                        process::exit(1);
                    }
                } else {
                    eprintln!("Error: -e necesita dos números: a inv_f");
                    process::exit(1);
                }
            }
            "-o" | "--output" => {
                i += 1;
                if i < args.len() {
                    output_base = args[i].clone();
                    i += 1;
                }
            }
            "-mstep" | "--max-step" => {
                i += 1;
                if i < args.len() {
                    if let Ok(step) = args[i].parse::<f64>() {
                        max_step = deg2rad(step);
                        i += 1;
                    }
                }
            }
            _ => {
                eprintln!("Argumento desconocido: {}", args[i]);
                i += 1;
            }
        }
    }

    // Validaciones
    if modo.is_empty() {
        eprintln!("Debe especificar -i, -d o -poly");
        process::exit(1);
    }

    let f = 1.0 / inv_f;
    let elli = Pelipsoide::new(f, semi_major);
    let plat = Platn6::new(elli.n);

    match modo.as_str() {
        "inverso" => {
            if p1.is_none() || p2.is_none() {
                eprintln!("El problema inverso requiere -P1 y -P2");
                process::exit(1);
            }
            let (lat1, lon1) = p1.unwrap();
            let (lat2, lon2) = p2.unwrap();
            let phi1 = deg2rad(lat1);
            let l1 = deg2rad(lon1);
            let phi2 = deg2rad(lat2);
            let l2 = deg2rad(lon2);

            let (alpha, dist, area, phi0, l0, pathpoints) = match tipo.as_str() {
                "align" => inv_align_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                "normal" => inv_normal_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                "geodesic" => inv_geodesic_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                "loxo" => inv_loxo_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                _ => inv_central_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
            };

            println!("\n========================================");
            println!("Acimut (deg): {:.10}", rad2deg(alpha));
            println!("Distancia (m): {:.4}", dist);
            println!("Area (m2): {:.0}", area);
            println!("Latitud Vértice phi0 (deg): {:.10}", rad2deg(phi0));
            println!("Longitud Vértice L0 (deg): {:.10}", rad2deg(l0));
            println!("========================================\n");

            if !output_base.is_empty() {
                let _ = guardar_kmz(&pathpoints, &format!("{}.kmz", output_base));
                let _ = guardar_shp(&pathpoints, &format!("{}.shp", output_base));
                let _ = guardar_shp_points(&pathpoints, &format!("{}.shp", output_base));
                let _ = guardar_csv(&pathpoints, &format!("{}.csv", output_base));
                println!("Archivos guardados con base: {}", output_base);
            }
        }
        "directo" => {
            if p1.is_none() || azimut.is_none() || distancia.is_none() {
                eprintln!("El problema directo requiere -P1, -az y -s");
                process::exit(1);
            }
            let (lat1, lon1) = p1.unwrap();
            let phi1 = deg2rad(lat1);
            let l1 = deg2rad(lon1);
            let alpha = deg2rad(azimut.unwrap());
            let dist = distancia.unwrap();

            let (phi2, l2) = match tipo.as_str() {
                "align" => direct_curva_align(&plat, &elli, phi1, l1, alpha, dist, max_step),
                "normal" => direct_curva_normal(&plat, &elli, phi1, l1, alpha, dist, max_step),
                "geodesic" => direct_curva_geodesic(&plat, &elli, phi1, l1, alpha, dist, max_step),
                "loxo" => direct_curva_loxo(&plat, &elli, phi1, l1, alpha, dist, max_step),
                _ => direct_curva_central(&plat, &elli, phi1, l1, alpha, dist, max_step),
            };
            let (alpha_calc, dist_calc, area, phi0, l0, pathpoints) = match tipo.as_str() {
                "align" => inv_align_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                "normal" => inv_normal_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                "geodesic" => inv_geodesic_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                "loxo" => inv_loxo_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
                _ => inv_central_dist_area(&plat, &elli, phi1, l1, phi2, l2, max_step),
            };

            println!("\n========================================");
            println!("Latitud P2 (deg): {:.10}", rad2deg(phi2));
            println!("Longitud P2 (deg): {:.10}", rad2deg(l2));
            println!("Acimut Calculado (deg): {:.10}", rad2deg(alpha_calc));
            println!("Distancia Calculada (m): {:.4}", dist_calc);
            println!("Area (m2): {:.0}", area);
            println!("Latitud Vértice phi0 (deg): {:.10}", rad2deg(phi0));
            println!("Longitud Vértice L0 (deg): {:.10}", rad2deg(l0));
            println!("========================================\n");

            if !output_base.is_empty() {
                // Para exportar la curva, calculamos el inverso hacia el punto encontrado
                let _ = guardar_kmz(&pathpoints, &format!("{}.kmz", output_base));
                let _ = guardar_shp(&pathpoints, &format!("{}.shp", output_base));
                let _ = guardar_shp_points(&pathpoints, &format!("{}.shp", output_base));
                let _ = guardar_csv(&pathpoints, &format!("{}.csv", output_base));
                println!("Archivos guardados con base: {}", output_base);
            }
        }
        "poly" => {
            let filename = poly_file.expect("Falta archivo para -poly");
            if let Ok(vertices) = leer_poligono(&filename) {
                let vertices_rad: Vec<(f64, f64)> = vertices.iter().map(|(lat, lon)| (deg2rad(*lat), deg2rad(*lon))).collect();

                let mut total_area = 0.0;
                let mut all_pathpoints: Vec<[f64; 5]> = Vec::new();

                println!("\n========================================");
                println!("Calculando polígono con {} vértices...", vertices_rad.len());

                for i in 0..vertices_rad.len() {
                    let (b1, l1) = vertices_rad[i];
                    let (b2, l2) = vertices_rad[(i + 1) % vertices_rad.len()];

                    let (alpha, dist, area, _, _, pathpoints) = match tipo.as_str() {
                        "align" => inv_align_dist_area(&plat, &elli, b1, l1, b2, l2, max_step),
                        "normal" => inv_normal_dist_area(&plat, &elli, b1, l1, b2, l2, max_step),
                        "geodesic" => inv_geodesic_dist_area(&plat, &elli, b1, l1, b2, l2, max_step),
                        "loxo" => inv_loxo_dist_area(&plat, &elli, b1, l1, b2, l2, max_step),
                        _ => inv_central_dist_area(&plat, &elli, b1, l1, b2, l2, max_step),
                    };

                    total_area += area;
                    if i == 0 {
                        all_pathpoints = pathpoints;
                    } else {
                        all_pathpoints.extend_from_slice(&pathpoints[1..]);
                    }

                    println!(
                        "Arista {} -> {}: Distancia = {:.4} m, Acimut = {:.8} deg, Área = {:.2} m²",
                        i + 1,
                        (i + 1) % vertices_rad.len() + 1,
                        dist,
                        rad2deg(alpha),
                        area
                    );
                }

                println!("----------------------------------------");
                println!("Superficie del polígono (m2): {:.2}", total_area.abs());
                println!("========================================\n");

                if !output_base.is_empty() {
                    let _ = guardar_kmz(&all_pathpoints, &format!("{}.kmz", output_base));
                    let _ = guardar_shp(&all_pathpoints, &format!("{}.shp", output_base));
                    let _ = guardar_csv(&all_pathpoints, &format!("{}.csv", output_base));
                    println!("Archivos guardados con base: {}", output_base);
                }
            } else {
                eprintln!("Error al leer el archivo de polígono: {}", filename);
                process::exit(1);
            }
        }
        _ => unreachable!(),
    }
}

fn print_help() {
    println!("Uso: curvas.exe [OPCIONES]");
    println!("Opciones:");
    println!("  -i, --inverso          Problema inverso (necesita -P1 y -P2)");
    println!("  -d, --directo          Problema directo (necesita -P1, -az y -s)");
    println!("  -poly, --poly-sup ARCHIVO  Área de polígono desde archivo CSV/TXT");
    println!("  -P1 LAT LON            Coordenadas del punto 1 (grados decimales)");
    println!("  -P2 LAT LON            Coordenadas del punto 2 (grados decimales)");
    println!("  -az, --azimuth GRADOS  Acimut inicial (para problema directo)");
    println!("  -s, --distance METROS  Distancia (para problema directo)");
    println!("  -t, --tipo TIPO        Tipo de curva o sección: align, central, normal, geodesic, loxo (defecto: central)");
    println!("  -o, --output BASE      Nombre base para archivos de salida (sin extensión)");
    println!("  -mstep, --max-step GRADOS  Paso máximo del integrador (defecto: 0.1°)");
}
