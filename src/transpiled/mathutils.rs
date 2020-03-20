extern "C" {
    #[no_mangle]
    fn cosf(_: f32) -> f32;
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn fvec_ones(s: *mut fvec_t);
    #[no_mangle]
    fn fvec_set_all(s: *mut fvec_t, val: smpl_t);
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    #[no_mangle]
    fn expf(_: f32) -> f32;
    #[no_mangle]
    fn logf(_: f32) -> f32;
    #[no_mangle]
    fn log10f(_: f32) -> f32;
    #[no_mangle]
    fn fabsf(_: f32) -> f32;
    #[no_mangle]
    fn powf(_: f32, _: f32) -> f32;
    #[no_mangle]
    fn floorf(_: f32) -> f32;
    #[no_mangle]
    fn strcmp(_: *const i8, _: *const i8) -> i32;
}
pub type smpl_t = f32;
pub type lsmp_t = f64;
pub type uint_t = u32;
pub type sint_t = i32;
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
/* * print format for sample in single precision */
/* * long sample format (64 bits or more) */
/* * print format for sample in double precision */
/* * unsigned integer */
/* * signed integer */
/* * character */
pub type char_t = i8;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
pub const aubio_win_parzen: aubio_window_type = 9;
pub const aubio_win_welch: aubio_window_type = 8;
pub const aubio_win_gaussian: aubio_window_type = 7;
pub const aubio_win_blackman_harris: aubio_window_type = 6;
pub const aubio_win_blackman: aubio_window_type = 5;
pub const aubio_win_hanningz: aubio_window_type = 4;
pub const aubio_win_hanning: aubio_window_type = 3;
pub const aubio_win_hamming: aubio_window_type = 2;
pub const aubio_win_rectangle: aubio_window_type = 1;
pub const aubio_win_ones: aubio_window_type = 0;
pub type aubio_window_type = u32;
pub const aubio_win_default: aubio_window_type = 4;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
pub type aubio_log_level = u32;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
#[no_mangle]
pub unsafe extern "C" fn new_aubio_window(window_type: *mut char_t,
                                          length: uint_t) -> *mut fvec_t {
    let win: *mut fvec_t = new_fvec(length);
    let mut err: uint_t = 0;
    if win.is_null() { return 0 as *mut fvec_t }
    err = fvec_set_window(win, window_type);
    if err != 0 as i32 as u32 {
        del_fvec(win);
        return 0 as *mut fvec_t
    }
    return win;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_set_window(win: *mut fvec_t,
                                         window_type: *mut char_t)
 -> uint_t {
    let w: *mut smpl_t = (*win).data;
    let mut i: uint_t = 0;
    let size: uint_t = (*win).length;
    let mut wintype: aubio_window_type = aubio_win_ones;
    if window_type.is_null() {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: window type can not be null.\n\x00" as
                      *const u8 as *const i8);
        return 1 as i32 as uint_t
    } else {
        if strcmp(window_type,
                  b"ones\x00" as *const u8 as *const i8) ==
               0 as i32 {
            wintype = aubio_win_ones
        } else if strcmp(window_type,
                         b"rectangle\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_rectangle
        } else if strcmp(window_type,
                         b"hamming\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_hamming
        } else if strcmp(window_type,
                         b"hanning\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_hanning
        } else if strcmp(window_type,
                         b"hanningz\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_hanningz
        } else if strcmp(window_type,
                         b"blackman\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_blackman
        } else if strcmp(window_type,
                         b"blackman_harris\x00" as *const u8 as
                             *const i8) == 0 as i32 {
            wintype = aubio_win_blackman_harris
        } else if strcmp(window_type,
                         b"gaussian\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_gaussian
        } else if strcmp(window_type,
                         b"welch\x00" as *const u8 as *const i8) ==
                      0 as i32 {
            wintype = aubio_win_welch
        } else if strcmp(window_type,
                         b"parzen\x00" as *const u8 as *const i8) ==
                      0 as i32 {
            wintype = aubio_win_parzen
        } else if strcmp(window_type,
                         b"default\x00" as *const u8 as *const i8)
                      == 0 as i32 {
            wintype = aubio_win_default
        } else {
            aubio_log(AUBIO_LOG_ERR as i32,
                      b"AUBIO ERROR: unknown window type `%s`.\n\x00" as
                          *const u8 as *const i8, window_type);
            return 1 as i32 as uint_t
        }
    }
    match wintype as u32 {
        0 => { fvec_ones(win); }
        1 => { fvec_set_all(win, 0.5f64 as smpl_t); }
        2 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (0.54f64 -
                         0.46f64 *
                             cosf((3.14159265358979323846264338327950288f64 *
                                       2.0f64 * i as f64 /
                                       size as f64) as
                                      f32) as f64) as
                        smpl_t;
                i = i.wrapping_add(1)
            }
        }
        3 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (0.5f64 -
                         0.5f64 *
                             cosf((3.14159265358979323846264338327950288f64 *
                                       2.0f64 * i as f64 /
                                       size as f64) as
                                      f32) as f64) as
                        smpl_t;
                i = i.wrapping_add(1)
            }
        }
        4 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (0.5f64 *
                         (1.0f64 -
                              cosf((3.14159265358979323846264338327950288f64 *
                                        2.0f64 * i as f64 /
                                        size as f64) as
                                       f32) as f64)) as
                        smpl_t;
                i = i.wrapping_add(1)
            }
        }
        5 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (0.42f64 -
                         0.50f64 *
                             cosf((3.14159265358979323846264338327950288f64 *
                                       2.0f64 * i as f64 /
                                       (size as f64 - 1.0f64)) as
                                      f32) as f64 +
                         0.08f64 *
                             cosf((2.0f64 *
                                       (3.14159265358979323846264338327950288f64
                                            * 2.0f64) * i as f64 /
                                       (size as f64 - 1.0f64)) as
                                      f32) as f64) as
                        smpl_t;
                i = i.wrapping_add(1)
            }
        }
        6 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (0.35875f64 -
                         0.48829f64 *
                             cosf((3.14159265358979323846264338327950288f64 *
                                       2.0f64 * i as f64 /
                                       (size as f64 - 1.0f64)) as
                                      f32) as f64 +
                         0.14128f64 *
                             cosf((2.0f64 *
                                       (3.14159265358979323846264338327950288f64
                                            * 2.0f64) * i as f64 /
                                       (size as f64 - 1.0f64)) as
                                      f32) as f64 -
                         0.01168f64 *
                             cosf((3.0f64 *
                                       (3.14159265358979323846264338327950288f64
                                            * 2.0f64) * i as f64 /
                                       (size as f64 - 1.0f64)) as
                                      f32) as f64) as
                        smpl_t;
                i = i.wrapping_add(1)
            }
        }
        7 => {
            let mut a: lsmp_t = 0.;
            let mut b: lsmp_t = 0.;
            let c: lsmp_t = 0.5f64;
            let mut n: uint_t = 0;
            n = 0 as i32 as uint_t;
            while n < size {
                a =
                    (n as f64 -
                         c *
                             size.wrapping_sub(1 as i32 as
                                                   u32) as
                                 f64) /
                        (c * c *
                             size.wrapping_sub(1 as i32 as
                                                   u32) as
                                 f64);
                b = -c * (a * a);
                *w.offset(n as isize) = expf(b as f32);
                n = n.wrapping_add(1)
            }
        }
        8 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (1.0f64 -
                         (2.0f64 * i as f64 -
                              size as f64) /
                             (size as f64 + 1.0f64) *
                             ((2.0f64 * i as f64 -
                                   size as f64) /
                                  (size as f64 + 1.0f64))) as
                        smpl_t;
                i = i.wrapping_add(1)
            }
        }
        9 => {
            i = 0 as i32 as uint_t;
            while i < size {
                *w.offset(i as isize) =
                    (1.0f64 -
                         fabsf((2.0f32 * i as f32 -
                                    size as f32) /
                                   (size as f32 + 1.0f32)) as
                             f64) as smpl_t;
                i = i.wrapping_add(1)
            }
        }
        _ => { }
    }
    return 0 as i32 as uint_t;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_unwrap2pi(phase: smpl_t) -> smpl_t {
    /* mod(phase+pi,-2pi)+pi */
    return (phase as f64 +
                3.14159265358979323846264338327950288f64 * 2.0f64 *
                    (1.0f64 +
                         floorf((-(phase as f64 +
                                       3.14159265358979323846264338327950288f64)
                                     /
                                     (3.14159265358979323846264338327950288f64
                                          * 2.0f64)) as f32) as
                             f64)) as smpl_t;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_mean(s: *mut fvec_t) -> smpl_t {
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        tmp += *(*s).data.offset(j as isize);
        j = j.wrapping_add(1)
    }
    return tmp / (*s).length as smpl_t;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_sum(s: *mut fvec_t) -> smpl_t {
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        tmp += *(*s).data.offset(j as isize);
        j = j.wrapping_add(1)
    }
    return tmp;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_max(s: *mut fvec_t) -> smpl_t {
    let mut j: uint_t = 0;
    let mut tmp: smpl_t = *(*s).data.offset(0 as i32 as isize);
    j = 1 as i32 as uint_t;
    while j < (*s).length {
        tmp =
            if tmp > *(*s).data.offset(j as isize) {
                tmp
            } else { *(*s).data.offset(j as isize) };
        j = j.wrapping_add(1)
    }
    return tmp;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_min(s: *mut fvec_t) -> smpl_t {
    let mut j: uint_t = 0;
    let mut tmp: smpl_t = *(*s).data.offset(0 as i32 as isize);
    j = 1 as i32 as uint_t;
    while j < (*s).length {
        tmp =
            if tmp < *(*s).data.offset(j as isize) {
                tmp
            } else { *(*s).data.offset(j as isize) };
        j = j.wrapping_add(1)
    }
    return tmp;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_min_elem(s: *mut fvec_t) -> uint_t {
    let mut j: uint_t = 0;
    let mut pos: uint_t = 0.0f64 as uint_t;
    let mut tmp: smpl_t = *(*s).data.offset(0 as i32 as isize);
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        pos = if tmp < *(*s).data.offset(j as isize) { pos } else { j };
        tmp =
            if tmp < *(*s).data.offset(j as isize) {
                tmp
            } else { *(*s).data.offset(j as isize) };
        j = j.wrapping_add(1)
    }
    return pos;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_max_elem(s: *mut fvec_t) -> uint_t {
    let mut j: uint_t = 0;
    let mut pos: uint_t = 0 as i32 as uint_t;
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        pos = if tmp > *(*s).data.offset(j as isize) { pos } else { j };
        tmp =
            if tmp > *(*s).data.offset(j as isize) {
                tmp
            } else { *(*s).data.offset(j as isize) };
        j = j.wrapping_add(1)
    }
    return pos;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_shift(s: *mut fvec_t) {
    let half: uint_t =
        (*s).length.wrapping_div(2 as i32 as u32);
    let mut start: uint_t = half;
    let mut j: uint_t = 0;
    // if length is odd, middle element is moved to the end
    if (2 as i32 as u32).wrapping_mul(half) < (*s).length {
        start = start.wrapping_add(1)
    }
    j = 0 as i32 as uint_t;
    while j < half {
        let t: smpl_t = *(*s).data.offset(j as isize);
        *(*s).data.offset(j as isize) =
            *(*s).data.offset(j.wrapping_add(start) as isize);
        *(*s).data.offset(j.wrapping_add(start) as isize) = t;
        j = j.wrapping_add(1)
    }
    if start != half {
        j = 0 as i32 as uint_t;
        while j < half {
            let t_0: smpl_t =
                *(*s).data.offset(j.wrapping_add(start).wrapping_sub(1 as
                                                                         i32
                                                                         as
                                                                         u32)
                                      as isize);
            *(*s).data.offset(j.wrapping_add(start).wrapping_sub(1 as
                                                                     i32
                                                                     as
                                                                     u32)
                                  as isize) =
                *(*s).data.offset(j.wrapping_add(start) as isize);
            *(*s).data.offset(j.wrapping_add(start) as isize) = t_0;
            j = j.wrapping_add(1)
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_ishift(s: *mut fvec_t) {
    let half: uint_t =
        (*s).length.wrapping_div(2 as i32 as u32);
    let mut start: uint_t = half;
    let mut j: uint_t = 0;
    // if length is odd, middle element is moved to the beginning
    if (2 as i32 as u32).wrapping_mul(half) < (*s).length {
        start = start.wrapping_add(1)
    }
    j = 0 as i32 as uint_t;
    while j < half {
        let t: smpl_t = *(*s).data.offset(j as isize);
        *(*s).data.offset(j as isize) =
            *(*s).data.offset(j.wrapping_add(start) as isize);
        *(*s).data.offset(j.wrapping_add(start) as isize) = t;
        j = j.wrapping_add(1)
    }
    if start != half {
        j = 0 as i32 as uint_t;
        while j < half {
            let t_0: smpl_t = *(*s).data.offset(half as isize);
            *(*s).data.offset(half as isize) = *(*s).data.offset(j as isize);
            *(*s).data.offset(j as isize) = t_0;
            j = j.wrapping_add(1)
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_push(in_0: *mut fvec_t,
                                   new_elem: smpl_t) {
    let mut i: uint_t = 0;
    i = 0 as i32 as uint_t;
    while i < (*in_0).length.wrapping_sub(1 as i32 as u32) {
        *(*in_0).data.offset(i as isize) =
            *(*in_0).data.offset(i.wrapping_add(1 as i32 as
                                                    u32) as isize);
        i = i.wrapping_add(1)
    }
    *(*in_0).data.offset((*in_0).length.wrapping_sub(1 as i32 as
                                                         u32) as
                             isize) = new_elem;
}
/* * clamp the values of a vector within the range [-abs(max), abs(max)]

  \param in vector to clamp
  \param absmax maximum value over which input vector elements should be clamped

*/
#[no_mangle]
pub unsafe extern "C" fn fvec_clamp(in_0: *mut fvec_t,
                                    absmax: smpl_t) {
    let mut i: uint_t = 0;
    i = 0 as i32 as uint_t;
    while i < (*in_0).length {
        if *(*in_0).data.offset(i as isize) >
               0 as i32 as f32 &&
               *(*in_0).data.offset(i as isize) > fabsf(absmax) {
            *(*in_0).data.offset(i as isize) = absmax
        } else if *(*in_0).data.offset(i as isize) <
                      0 as i32 as f32 &&
                      *(*in_0).data.offset(i as isize) < -fabsf(absmax) {
            *(*in_0).data.offset(i as isize) = -absmax
        }
        i = i.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_level_lin(f: *const fvec_t) -> smpl_t {
    let mut energy: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*f).length {
        energy +=
            *(*f).data.offset(j as isize) * *(*f).data.offset(j as isize);
        j = j.wrapping_add(1)
    }
    return energy / (*f).length as f32;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_local_hfc(v: *mut fvec_t) -> smpl_t {
    let mut hfc: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*v).length {
        hfc +=
            j.wrapping_add(1 as i32 as u32) as f32
                * *(*v).data.offset(j as isize);
        j = j.wrapping_add(1)
    }
    return hfc;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_min_removal(v: *mut fvec_t) {
    let v_min: smpl_t = fvec_min(v);
    fvec_add(v, -v_min);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_alpha_norm(o: *mut fvec_t,
                                         alpha: smpl_t) -> smpl_t {
    let mut j: uint_t = 0;
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*o).length {
        tmp += powf(fabsf(*(*o).data.offset(j as isize)), alpha);
        j = j.wrapping_add(1)
    }
    return powf(tmp / (*o).length as f32,
                (1.0f64 / alpha as f64) as f32);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_alpha_normalise(o: *mut fvec_t,
                                              alpha: smpl_t) {
    let mut j: uint_t = 0;
    let norm: smpl_t = fvec_alpha_norm(o, alpha);
    j = 0 as i32 as uint_t;
    while j < (*o).length {
        let ref mut fresh0 = *(*o).data.offset(j as isize);
        *fresh0 /= norm;
        j = j.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_add(o: *mut fvec_t, val: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*o).length {
        let ref mut fresh1 = *(*o).data.offset(j as isize);
        *fresh1 += val;
        j = j.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_mul(o: *mut fvec_t, val: smpl_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*o).length {
        let ref mut fresh2 = *(*o).data.offset(j as isize);
        *fresh2 *= val;
        j = j.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_adapt_thres(vec: *mut fvec_t,
                                          tmp: *mut fvec_t,
                                          post: uint_t, pre: uint_t) {
    let length: uint_t = (*vec).length;
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < length {
        let ref mut fresh3 = *(*vec).data.offset(j as isize);
        *fresh3 -= fvec_moving_thres(vec, tmp, post, pre, j);
        j = j.wrapping_add(1)
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_moving_thres(vec: *mut fvec_t,
                                           tmpvec: *mut fvec_t,
                                           post: uint_t, pre: uint_t,
                                           pos: uint_t) -> smpl_t {
    let mut k: uint_t = 0;
    let medar: *mut smpl_t = (*tmpvec).data;
    let win_length: uint_t =
        post.wrapping_add(pre).wrapping_add(1 as i32 as u32);
    let length: uint_t = (*vec).length;
    /* post part of the buffer does not exist */
    if pos < post.wrapping_add(1 as i32 as u32) {
        k = 0 as i32 as uint_t; /* 0-padding at the beginning */
        while k <
                  post.wrapping_add(1 as i32 as
                                        u32).wrapping_sub(pos) {
            *medar.offset(k as isize) = 0.0f64 as smpl_t;
            k = k.wrapping_add(1)
        }
        k =
            post.wrapping_add(1 as i32 as
                                  u32).wrapping_sub(pos);
        while k < win_length {
            *medar.offset(k as isize) =
                *(*vec).data.offset(k.wrapping_add(pos).wrapping_sub(post) as
                                        isize);
            k = k.wrapping_add(1)
        }
        /* the buffer is fully defined */
    } else if pos.wrapping_add(pre) < length {
        k = 0 as i32 as uint_t;
        while k < win_length {
            *medar.offset(k as isize) =
                *(*vec).data.offset(k.wrapping_add(pos).wrapping_sub(post) as
                                        isize);
            k = k.wrapping_add(1)
        }
        /* pre part of the buffer does not exist */
    } else {
        k = 0 as i32 as uint_t;
        while k < length.wrapping_sub(pos).wrapping_add(post) {
            *medar.offset(k as isize) =
                *(*vec).data.offset(k.wrapping_add(pos).wrapping_sub(post) as
                                        isize);
            k = k.wrapping_add(1)
        }
        k = length.wrapping_sub(pos).wrapping_add(post);
        while k < win_length {
            *medar.offset(k as isize) = 0.0f64 as smpl_t;
            k = k.wrapping_add(1)
        }
        /* 0-padding at the end */
    }
    return fvec_median(tmpvec);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_median(input: *mut fvec_t) -> smpl_t {
    let n: uint_t = (*input).length;
    let arr: *mut smpl_t = (*input).data;
    let mut low: uint_t = 0;
    let mut high: uint_t = 0;
    let mut median: uint_t = 0;
    let mut middle: uint_t = 0;
    let mut ll: uint_t = 0;
    let mut hh: uint_t = 0;
    low = 0 as i32 as uint_t;
    high = n.wrapping_sub(1 as i32 as u32);
    median =
        low.wrapping_add(high).wrapping_div(2 as i32 as u32);
    loop  {
        if high <= low {
            /* One element only */
            return *arr.offset(median as isize)
        }
        if high == low.wrapping_add(1 as i32 as u32) {
            /* Two elements only */
            if *arr.offset(low as isize) > *arr.offset(high as isize) {
                let t: smpl_t = *arr.offset(low as isize);
                *arr.offset(low as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t
            }
            return *arr.offset(median as isize)
        }
        /* Find median of low, middle and high items; swap into position low */
        middle =
            low.wrapping_add(high).wrapping_div(2 as i32 as
                                                    u32);
        if *arr.offset(middle as isize) > *arr.offset(high as isize) {
            let t_0: smpl_t = *arr.offset(middle as isize);
            *arr.offset(middle as isize) = *arr.offset(high as isize);
            *arr.offset(high as isize) = t_0
        }
        if *arr.offset(low as isize) > *arr.offset(high as isize) {
            let t_1: smpl_t = *arr.offset(low as isize);
            *arr.offset(low as isize) = *arr.offset(high as isize);
            *arr.offset(high as isize) = t_1
        }
        if *arr.offset(middle as isize) > *arr.offset(low as isize) {
            let t_2: smpl_t = *arr.offset(middle as isize);
            *arr.offset(middle as isize) = *arr.offset(low as isize);
            *arr.offset(low as isize) = t_2
        }
        /* Swap low item (now in position middle) into position (low+1) */
        let t_3: smpl_t = *arr.offset(middle as isize);
        *arr.offset(middle as isize) =
            *arr.offset(low.wrapping_add(1 as i32 as u32) as
                            isize);
        *arr.offset(low.wrapping_add(1 as i32 as u32) as
                        isize) = t_3;
        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low.wrapping_add(1 as i32 as u32);
        hh = high;
        loop  {
            loop  {
                ll = ll.wrapping_add(1);
                if !(*arr.offset(low as isize) > *arr.offset(ll as isize)) {
                    break ;
                }
            }
            loop  {
                hh = hh.wrapping_sub(1);
                if !(*arr.offset(hh as isize) > *arr.offset(low as isize)) {
                    break ;
                }
            }
            if hh < ll { break ; }
            let t_4: smpl_t = *arr.offset(ll as isize);
            *arr.offset(ll as isize) = *arr.offset(hh as isize);
            *arr.offset(hh as isize) = t_4
        }
        /* Swap middle item (in position low) back into correct position */
        let t_5: smpl_t = *arr.offset(low as isize);
        *arr.offset(low as isize) = *arr.offset(hh as isize);
        *arr.offset(hh as isize) = t_5;
        /* Re-set active partition */
        if hh <= median { low = ll } // avoid nans and infs
        if hh >= median {
            high = hh.wrapping_sub(1 as i32 as u32)
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn fvec_quadratic_peak_pos(x: *const fvec_t,
                                                 pos: uint_t) -> smpl_t {
    let mut s0: smpl_t = 0.;
    let mut s1: smpl_t = 0.;
    let mut s2: smpl_t = 0.;
    let mut x0: uint_t = 0;
    let mut x2: uint_t = 0;
    let half: smpl_t = 0.5f64 as smpl_t;
    let two: smpl_t = 2.0f64 as smpl_t;
    if pos == 0 as i32 as u32 ||
           pos == (*x).length.wrapping_sub(1 as i32 as u32) {
        return pos as smpl_t
    }
    x0 =
        if pos < 1 as i32 as u32 {
            pos
        } else { pos.wrapping_sub(1 as i32 as u32) };
    x2 =
        if pos.wrapping_add(1 as i32 as u32) < (*x).length {
            pos.wrapping_add(1 as i32 as u32)
        } else { pos };
    if x0 == pos {
        return if *(*x).data.offset(pos as isize) <=
                      *(*x).data.offset(x2 as isize) {
                   pos
               } else { x2 } as smpl_t
    }
    if x2 == pos {
        return if *(*x).data.offset(pos as isize) <=
                      *(*x).data.offset(x0 as isize) {
                   pos
               } else { x0 } as smpl_t
    }
    s0 = *(*x).data.offset(x0 as isize);
    s1 = *(*x).data.offset(pos as isize);
    s2 = *(*x).data.offset(x2 as isize);
    return pos as f32 + half * (s0 - s2) / (s0 - two * s1 + s2);
}
#[no_mangle]
pub unsafe extern "C" fn fvec_quadratic_peak_mag(x: *mut fvec_t,
                                                 pos: smpl_t) -> smpl_t {
    let mut x0: smpl_t = 0.;
    let mut x1: smpl_t = 0.;
    let mut x2: smpl_t = 0.;
    let index: uint_t =
        ((pos as f64 - 0.5f64) as
             uint_t).wrapping_add(1 as i32 as u32);
    if pos >= (*x).length as f32 || (pos as f64) < 0.0f64
       {
        return 0.0f64 as smpl_t
    }
    if index as smpl_t == pos { return *(*x).data.offset(index as isize) }
    x0 =
        *(*x).data.offset(index.wrapping_sub(1 as i32 as u32)
                              as isize);
    x1 = *(*x).data.offset(index as isize);
    x2 =
        *(*x).data.offset(index.wrapping_add(1 as i32 as u32)
                              as isize);
    return (x1 as f64 -
                0.25f64 * (x0 - x2) as f64 *
                    (pos - index as f32) as f64) as
               smpl_t;
}
#[no_mangle]
pub unsafe extern "C" fn fvec_peakpick(onset: *const fvec_t,
                                       pos: uint_t) -> uint_t {
    let mut tmp: uint_t = 0 as i32 as uint_t;
    tmp =
        (*(*onset).data.offset(pos as isize) >
             *(*onset).data.offset(pos.wrapping_sub(1 as i32 as
                                                        u32) as
                                       isize) &&
             *(*onset).data.offset(pos as isize) >
                 *(*onset).data.offset(pos.wrapping_add(1 as i32 as
                                                            u32) as
                                           isize) &&
             *(*onset).data.offset(pos as isize) as f64 > 0.0f64)
            as i32 as uint_t;
    return tmp;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_quadfrac(s0: smpl_t, s1: smpl_t,
                                        s2: smpl_t, pf: smpl_t)
 -> smpl_t {
    let tmp: smpl_t =
        (s0 as f64 +
             pf as f64 / 2.0f64 *
                 (pf as f64 *
                      (s0 as f64 - 2.0f64 * s1 as f64 +
                           s2 as f64) -
                      3.0f64 * s0 as f64 +
                      4.0f64 * s1 as f64 - s2 as f64))
            as smpl_t;
    return tmp;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_freqtomidi(freq: smpl_t) -> smpl_t {
    let mut midi: smpl_t = 0.;
    if (freq as f64) < 2.0f64 ||
           freq as f64 > 100000.0f64 {
        return 0.0f64 as smpl_t
    }
    /* log(freq/A-2)/log(2) */
    midi = (freq as f64 / 6.875f64) as smpl_t; // avoid infs
    midi = (logf(midi) as f64 / 0.6931471805599453f64) as smpl_t;
    midi *= 12 as i32 as f32;
    midi -= 3 as i32 as f32;
    return midi;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_miditofreq(midi: smpl_t) -> smpl_t {
    let mut freq: smpl_t = 0.;
    if midi as f64 > 140.0f64 { return 0.0f64 as smpl_t }
    freq = ((midi as f64 + 3.0f64) / 12.0f64) as smpl_t;
    freq =
        expf((freq as f64 * 0.6931471805599453f64) as
                 f32);
    freq = (freq as f64 * 6.875f64) as smpl_t;
    return freq;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_bintofreq(bin: smpl_t,
                                         samplerate: smpl_t,
                                         fftsize: smpl_t) -> smpl_t {
    let freq: smpl_t = samplerate / fftsize;
    return freq *
               (if bin > 0 as i32 as f32 {
                    bin
                } else { 0 as i32 as f32 });
}
#[no_mangle]
pub unsafe extern "C" fn aubio_bintomidi(bin: smpl_t,
                                         samplerate: smpl_t,
                                         fftsize: smpl_t) -> smpl_t {
    let midi: smpl_t = aubio_bintofreq(bin, samplerate, fftsize);
    return aubio_freqtomidi(midi);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_freqtobin(freq: smpl_t,
                                         samplerate: smpl_t,
                                         fftsize: smpl_t) -> smpl_t {
    let bin: smpl_t = fftsize / samplerate;
    return (if freq > 0 as i32 as f32 {
                freq
            } else { 0 as i32 as f32 }) * bin;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_miditobin(midi: smpl_t,
                                         samplerate: smpl_t,
                                         fftsize: smpl_t) -> smpl_t {
    let freq: smpl_t = aubio_miditofreq(midi);
    return aubio_freqtobin(freq, samplerate, fftsize);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_is_power_of_two(a: uint_t) -> uint_t {
    if a & a.wrapping_sub(1 as i32 as u32) ==
           0 as i32 as u32 {
        return 1 as i32 as uint_t
    } else { return 0 as i32 as uint_t };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_next_power_of_two(a: uint_t) -> uint_t {
    let mut i: uint_t = 1 as i32 as uint_t;
    while i < a { i <<= 1 as i32 }
    return i;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_power_of_two_order(a: uint_t) -> uint_t {
    let mut order: i32 = 0 as i32;
    let mut temp: i32 = aubio_next_power_of_two(a) as i32;
    loop  {
        temp >>= 1 as i32;
        if !(temp != 0) { break ; }
        order += 1
    }
    return order as uint_t;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_db_spl(o: *const fvec_t) -> smpl_t {
    return (10.0f64 * log10f(aubio_level_lin(o)) as f64) as smpl_t;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_silence_detection(o: *const fvec_t,
                                                 threshold: smpl_t)
 -> uint_t {
    return (aubio_db_spl(o) < threshold) as i32 as uint_t;
}
/* * compute sound level on a linear scale

  This gives the average of the square amplitudes.

  \param v vector to compute level from

  \return level of v

*/
/* * compute sound pressure level (SPL) in dB

  This quantity is often wrongly called 'loudness'.

  This gives ten times the log10 of the average of the square amplitudes.

  \param v vector to compute dB SPL from

  \return level of v in dB SPL

*/
/* * check if buffer level in dB SPL is under a given threshold

  \param v vector to get level from
  \param threshold threshold in dB SPL

  \return 0 if level is under the given threshold, 1 otherwise

*/
/* * get buffer level if level >= threshold, 1. otherwise

  \param v vector to get level from
  \param threshold threshold in dB SPL

  \return level in dB SPL if level >= threshold, 1. otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_level_detection(o: *const fvec_t,
                                               threshold: smpl_t)
 -> smpl_t {
    let db_spl: smpl_t = aubio_db_spl(o);
    if db_spl < threshold { return 1.0f64 as smpl_t } else { return db_spl };
}
/* * zero-crossing rate (ZCR)

  The zero-crossing rate is the number of times a signal changes sign,
  divided by the length of this signal.

  \param v vector to compute ZCR from

  \return zero-crossing rate of v

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_zero_crossing_rate(input: *mut fvec_t)
 -> smpl_t {
    let mut j: uint_t = 0;
    let mut zcr: uint_t = 0 as i32 as uint_t;
    j = 1 as i32 as uint_t;
    while j < (*input).length {
        // previous was strictly negative
        if (*(*input).data.offset(j.wrapping_sub(1 as i32 as
                                                     u32) as isize)
                as f64) < 0.0f64 {
            // current is positive or null
            if *(*input).data.offset(j as isize) as f64 >= 0.0f64 {
                zcr =
                    (zcr as
                         u32).wrapping_add(1 as i32 as
                                                        u32) as
                        uint_t as uint_t
            }
            // previous was positive or null
        } else if (*(*input).data.offset(j as isize) as f64) <
                      0.0f64 {
            zcr =
                (zcr as
                     u32).wrapping_add(1 as i32 as
                                                    u32) as uint_t as
                    uint_t
        }
        j = j.wrapping_add(1)
    }
    return zcr as f32 / (*input).length as smpl_t;
}
// current is strictly negative
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Various math functions

  \example test-mathutils.c
  \example test-mathutils-window.c

 */
/* * compute the mean of a vector

  \param s vector to compute mean from
  \return the mean of `v`

*/
/* * find the max of a vector

  \param s vector to get the max from

  \return the value of the minimum of v

*/
/* * find the min of a vector

  \param s vector to get the min from

  \return the value of the maximum of v

*/
/* * find the index of the min of a vector

  \param s vector to get the index from

  \return the index of the minimum element of v

*/
/* * find the index of the max of a vector

  \param s vector to get the index from

  \return the index of the maximum element of v

*/
/* * swap the left and right halves of a vector

  This function swaps the left part of the signal with the right part of the
signal. Therefore

  \f$ a[0], a[1], ..., a[\frac{N}{2}], a[\frac{N}{2}+1], ..., a[N-1], a[N] \f$

  becomes

  \f$ a[\frac{N}{2}+1], ..., a[N-1], a[N], a[0], a[1], ..., a[\frac{N}{2}] \f$

  This operation, known as 'fftshift' in the Matlab Signal Processing Toolbox,
can be used before computing the FFT to simplify the phase relationship of the
resulting spectrum. See Amalia de Götzen's paper referred to above.

*/
/* * swap the left and right halves of a vector

  This function swaps the left part of the signal with the right part of the
signal. Therefore

  \f$ a[0], a[1], ..., a[\frac{N}{2}], a[\frac{N}{2}+1], ..., a[N-1], a[N] \f$

  becomes

  \f$ a[\frac{N}{2}+1], ..., a[N-1], a[N], a[0], a[1], ..., a[\frac{N}{2}] \f$

  This operation, known as 'ifftshift' in the Matlab Signal Processing Toolbox,
can be used after computing the inverse FFT to simplify the phase relationship
of the resulting spectrum. See Amalia de Götzen's paper referred to above.

*/
/* * push a new element to the end of a vector, erasing the first element and
 * sliding all others

  \param in vector to push to
  \param new_elem new_element to add at the end of the vector

  In numpy words, this is equivalent to: in = np.concatenate([in, [new_elem]])[1:]

*/
/* * compute the sum of all elements of a vector

  \param v vector to compute the sum of

  \return the sum of v

*/
/* * compute the High Frequency Content of a vector

  The High Frequency Content is defined as \f$ \sum_0^{N-1} (k+1) v[k] \f$.

  \param v vector to get the energy from

  \return the HFC of v

*/
/* * computes the p-norm of a vector

  Computes the p-norm of a vector for \f$ p = \alpha \f$

  \f$ L^p = ||x||_p = (|x_1|^p + |x_2|^p + ... + |x_n|^p ) ^ \frac{1}{p} \f$

  If p = 1, the result is the Manhattan distance.

  If p = 2, the result is the Euclidean distance.

  As p tends towards large values, \f$ L^p \f$ tends towards the maximum of the
input vector.

  References:

    - <a href="http://en.wikipedia.org/wiki/Lp_space">\f$L^p\f$ space</a> on
  Wikipedia

  \param v vector to compute norm from
  \param p order of the computed norm

  \return the p-norm of v

*/
/* *  alpha normalisation

  This function divides all elements of a vector by the p-norm as computed by
fvec_alpha_norm().

  \param v vector to compute norm from
  \param p order of the computed norm

*/
/* * add a constant to each elements of a vector

  \param v vector to add constant to
  \param c constant to add to v

*/
/* * multiply each elements of a vector by a scalar

  \param v vector to add constant to
  \param s constant to scale v with

*/
/* * remove the minimum value of the vector to each elements

  \param v vector to remove minimum from

*/
/* * compute moving median threshold of a vector

  This function computes the moving median threshold value of at the given
position of a vector, taking the median among post elements before and up to
pre elements after pos.

  \param v input vector
  \param tmp temporary vector of length post+1+pre
  \param post length of causal part to take before pos
  \param pre length of anti-causal part to take after pos
  \param pos index to compute threshold for

  \return moving median threshold value

*/
/* * apply adaptive threshold to a vector

  For each points at position p of an input vector, this function remove the
moving median threshold computed at p.

  \param v input vector
  \param tmp temporary vector of length post+1+pre
  \param post length of causal part to take before pos
  \param pre length of anti-causal part to take after pos

*/
/* * returns the median of a vector

  The QuickSelect routine is based on the algorithm described in "Numerical
recipes in C", Second Edition, Cambridge University Press, 1992, Section 8.5,
ISBN 0-521-43108-5

  This implementation of the QuickSelect routine is based on Nicolas
Devillard's implementation, available at http://ndevilla.free.fr/median/median/
and in the Public Domain.

  \param v vector to get median from

  \return the median of v

*/
/* * finds exact peak index by quadratic interpolation

  See [Quadratic Interpolation of Spectral
  Peaks](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Peak_Interpolation.html),
  by Julius O. Smith III

  \f$ p_{frac} = \frac{1}{2} \frac {x[p-1] - x[p+1]} {x[p-1] - 2 x[p] + x[p+1]} \in [ -.5, .5] \f$

  \param x vector to get the interpolated peak position from
  \param p index of the peak in vector `x`
  \return \f$ p + p_{frac} \f$ exact peak position of interpolated maximum or minimum

*/
/* * finds magnitude of peak by quadratic interpolation

  See [Quadratic Interpolation of Spectral
  Peaks](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Peak_Interpolation.html),
  by Julius O. Smith III

  \param x vector to get the magnitude of the interpolated peak position from
  \param p index of the peak in vector `x`
  \return magnitude of interpolated peak

*/
/* * Quadratic interpolation using Lagrange polynomial.

  Inspired from ``Comparison of interpolation algorithms in real-time sound
processing'', Vladimir Arnost,

  \param s0,s1,s2 are 3 consecutive samples of a curve
  \param pf is the floating point index [0;2]

  \return \f$ s0 + (pf/2.)*((pf-3.)*s0-2.*(pf-2.)*s1+(pf-1.)*s2); \f$

*/
/* * return 1 if v[p] is a peak and positive, 0 otherwise

  This function returns 1 if a peak is found at index p in the vector v. The
peak is defined as follows:

  - v[p] is positive
  - v[p-1] < v[p]
  - v[p] > v[p+1]

  \param v input vector
  \param p position of supposed for peak

  \return 1 if a peak is found, 0 otherwise

*/
/* * return 1 if a is a power of 2, 0 otherwise */
/* * return the next power of power of 2 greater than a */
/* * return the log2 factor of the given power of 2 value a */
/* * compute normalised autocorrelation function

  \param input vector to compute autocorrelation from
  \param output vector to store autocorrelation function to

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_autocorr(input: *const fvec_t,
                                        output: *mut fvec_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    let length: uint_t = (*input).length;
    let mut data: *mut smpl_t = 0 as *mut smpl_t;
    let mut acf: *mut smpl_t = 0 as *mut smpl_t;
    let mut tmp: smpl_t = 0 as i32 as smpl_t;
    data = (*input).data;
    acf = (*output).data;
    i = 0 as i32 as uint_t;
    while i < length {
        tmp = 0.0f64 as smpl_t;
        j = i;
        while j < length {
            tmp +=
                *data.offset(j.wrapping_sub(i) as isize) *
                    *data.offset(j as isize);
            j = j.wrapping_add(1)
        }
        *acf.offset(i as isize) = tmp / length.wrapping_sub(i) as smpl_t;
        i = i.wrapping_add(1)
    };
}
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * @file
 *  various functions useful in audio signal processing
 */
/* * create window

  \param window_type type of the window to create
  \param size length of the window to create (see fvec_set_window())

*/
/* * set elements of a vector to window coefficients

  \param window exsting ::fvec_t to use
  \param window_type type of the window to create

  List of available window types: "rectangle", "hamming", "hanning",
  "hanningz", "blackman", "blackman_harris", "gaussian", "welch", "parzen",
  "default".

  "default" is equivalent to "hanningz".

  References:

    - <a href="http://en.wikipedia.org/wiki/Window_function">Window
function</a> on Wikipedia
    - Amalia de Götzen, Nicolas Bernardini, and Daniel Arfib. Traditional (?)
implementations of a phase vocoder: the tricks of the trade. In Proceedings of
the International Conference on Digital Audio Effects (DAFx-00), pages 37–44,
Uni- versity of Verona, Italy, 2000.
  (<a href="http://www.cs.princeton.edu/courses/archive/spr09/cos325/Bernardini.pdf">
  pdf</a>)

 */
/* * compute the principal argument

  This function maps the input phase to its corresponding value wrapped in the
range \f$ [-\pi, \pi] \f$.

  \param phase unwrapped phase to map to the unit circle

  \return equivalent phase wrapped to the unit circle

*/
/* * convert frequency bin to midi value */
/* * convert midi value to frequency bin */
/* * convert frequency bin to frequency (Hz) */
/* * convert frequency (Hz) to frequency bin */
/* * convert frequency (Hz) to mel

  \param freq input frequency, in Hz

  \return output mel

  Converts a scalar from the frequency domain to the mel scale using Slaney
  Auditory Toolbox's implementation:

  If \f$ f < 1000 \f$, \f$ m = 3 f / 200 \f$.

  If \f$ f >= 1000 \f$, \f$ m = 1000 + 27 \frac{{ln}(f) - ln(1000))}
  {{ln}(6400) - ln(1000)}
  \f$

  See also
  --------

  aubio_meltohz(), aubio_hztomel_htk().

*/
/* * convert mel to frequency (Hz)

  \param mel input mel

  \return output frequency, in Hz

  Converts a scalar from the mel scale to the frequency domain using Slaney
  Auditory Toolbox's implementation:

  If \f$ f < 1000 \f$, \f$ f = 200 m/3 \f$.

  If \f$ f \geq 1000 \f$, \f$ f = 1000 + \left(\frac{6400}{1000}\right)
  ^{\frac{m - 1000}{27}} \f$

  See also
  --------

  aubio_hztomel(), aubio_meltohz_htk().

  References
  ----------

  Malcolm Slaney, *Auditory Toolbox Version 2, Technical Report #1998-010*
  https://engineering.purdue.edu/~malcolm/interval/1998-010/

*/
/* * convert frequency (Hz) to mel

  \param freq input frequency, in Hz

  \return output mel

  Converts a scalar from the frequency domain to the mel scale, using the
  equation defined by O'Shaughnessy, as implemented in the HTK speech
  recognition toolkit:

  \f$ m = 1127 + ln(1 + \frac{f}{700}) \f$

  See also
  --------

  aubio_meltohz_htk(), aubio_hztomel().

  References
  ----------

  Douglas O'Shaughnessy (1987). *Speech communication: human and machine*.
  Addison-Wesley. p. 150. ISBN 978-0-201-16520-3.

  HTK Speech Recognition Toolkit: http://htk.eng.cam.ac.uk/

 */
/* * convert mel to frequency (Hz)

  \param mel input mel

  \return output frequency, in Hz

  Converts a scalar from the mel scale to the frequency domain, using the
  equation defined by O'Shaughnessy, as implemented in the HTK speech
  recognition toolkit:

  \f$ f = 700 * {e}^\left(\frac{f}{1127} - 1\right) \f$

  See also
  --------

  aubio_hztomel_htk(), aubio_meltohz().

*/
/* * convert frequency (Hz) to midi value (0-128) */
/* * convert midi value (0-128) to frequency (Hz) */
/* * clean up cached memory at the end of program

  This function should be used at the end of programs to purge all cached
  memory. So far it is only useful to clean FFTW's cache.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_cleanup() { }
