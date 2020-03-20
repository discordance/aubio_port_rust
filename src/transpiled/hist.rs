use crate::transpiled::scale::aubio_scale_t;

extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn floorf(_: f32) -> f32;
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_zeros(s: *mut fvec_t);
    /* * find the max of a vector

  \param s vector to get the max from

  \return the value of the minimum of v

*/
    #[no_mangle]
    fn fvec_max(s: *mut fvec_t) -> smpl_t;
    /* * find the min of a vector

  \param s vector to get the min from

  \return the value of the maximum of v

*/
    #[no_mangle]
    fn fvec_min(s: *mut fvec_t) -> smpl_t;
    /* * create a scale object

  \param flow lower value of output function
  \param fhig higher value of output function
  \param ilow lower value of input function
  \param ihig higher value of output function

*/
    #[no_mangle]
    fn new_aubio_scale(flow: smpl_t, fhig: smpl_t, ilow: smpl_t, ihig: smpl_t)
     -> *mut aubio_scale_t;
    /* * delete a scale object

  \param s scale object as returned by new_aubio_scale

*/
    #[no_mangle]
    fn del_aubio_scale(s: *mut aubio_scale_t);
    /* * scale input vector

  \param s scale object as returned by new_aubio_scale
  \param input vector to scale

*/
    #[no_mangle]
    fn aubio_scale_do(s: *mut aubio_scale_t, input: *mut fvec_t);
    /* * modify scale parameters after object creation

  \param s scale object as returned by new_aubio_scale
  \param olow lower value of output function
  \param ohig higher value of output function
  \param ilow lower value of input function
  \param ihig higher value of output function

*/
    #[no_mangle]
    fn aubio_scale_set_limits(s: *mut aubio_scale_t, ilow: smpl_t,
                              ihig: smpl_t, olow: smpl_t, ohig: smpl_t)
     -> uint_t;
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
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
pub type smpl_t = f32;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
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

  Vector of real-valued data

  This file specifies the ::fvec_t buffer type, which is used throughout aubio
  to store vector of real-valued ::smpl_t.

  \example test-fvec.c

*/
/* * Buffer for real data

  Vector of real-valued data

  ::fvec_t is is the structure used to store vector of real-valued data, ::smpl_t .

  \code

  uint_t buffer_size = 1024;

  // create a vector of 512 values
  fvec_t * input = new_fvec (buffer_size);

  // set some values of the vector
  input->data[23] = 2.;
  // ..

  // compute the mean of the vector
  mean = fvec_mean(a_vector);

  // destroy the vector
  del_fvec(a_vector);

  \endcode

  See `examples/` and `tests/src` directories for more examples.

 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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

 Vector scaling function

 This object, inspired from the scale object in FTS, the jMax engine, scales
 the values of a vector according to an affine function defined as follow:

 \f$ y = (x - ilow)*(ohig-olow)/(ihig-ilow) + olow \f$

*/
/* * scale object */
// pub type aubio_scale_t = _aubio_scale_t;
/*
  Copyright (C) 2003-2009 Paul Brossier <piem@aubio.org>

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
//fvec_min fvec_max
/* *******
 * Object Structure
 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_hist_t {
    pub hist: *mut fvec_t,
    pub nelems: uint_t,
    pub cent: *mut fvec_t,
    pub scaler: *mut aubio_scale_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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
 *
 * Histogram function
 *
 * Big hacks to implement an histogram
 */
/* * histogram object */
pub type aubio_hist_t = _aubio_hist_t;
/* * histogram creation

  \param flow minimum input
  \param fhig maximum input
  \param nelems number of histogram columns

*/
/* *
 * Object creation/deletion calls
 */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_hist(flow: smpl_t, fhig: smpl_t,
                                        nelems: uint_t)
 -> *mut aubio_hist_t {
    let mut s: *mut aubio_hist_t =
        calloc(::std::mem::size_of::<aubio_hist_t>() as u64,
               1 as i32 as u64) as *mut aubio_hist_t;
    let step: smpl_t = (fhig - flow) / nelems as smpl_t;
    let mut accum: smpl_t = step;
    let mut i: uint_t = 0;
    if nelems as sint_t <= 0 as i32 {
        free(s as *mut core::ffi::c_void);
        return 0 as *mut aubio_hist_t
    }
    (*s).nelems = nelems;
    (*s).hist = new_fvec(nelems);
    (*s).cent = new_fvec(nelems);
    /* use scale to map flow/fhig -> 0/nelems */
    (*s).scaler =
        new_aubio_scale(flow, fhig, 0 as i32 as smpl_t,
                        nelems as smpl_t);
    /* calculate centers now once */
    *(*(*s).cent).data.offset(0 as i32 as isize) =
        (flow as f64 + 0.5f64 * step as f64) as smpl_t;
    i = 1 as i32 as uint_t;
    while i < (*s).nelems {
        *(*(*s).cent).data.offset(i as isize) =
            *(*(*s).cent).data.offset(0 as i32 as isize) + accum;
        i = i.wrapping_add(1);
        accum += step
    }
    return s;
}
/* * histogram deletion */
#[no_mangle]
pub unsafe extern "C" fn del_aubio_hist(s: *mut aubio_hist_t) {
    del_fvec((*s).hist);
    del_fvec((*s).cent);
    del_aubio_scale((*s).scaler);
    free(s as *mut core::ffi::c_void);
}
/* * compute the histogram */
/* **
 * do it
 */
#[no_mangle]
pub unsafe extern "C" fn aubio_hist_do(s: *mut aubio_hist_t,
                                       input: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut tmp: sint_t = 0 as i32;
    aubio_scale_do((*s).scaler, input);
    /* reset data */
    fvec_zeros((*s).hist);
    /* run accum */
    j = 0 as i32 as uint_t;
    while j < (*input).length {
        tmp = floorf(*(*input).data.offset(j as isize)) as sint_t;
        if tmp >= 0 as i32 && tmp < (*s).nelems as sint_t {
            let ref mut fresh0 = *(*(*s).hist).data.offset(tmp as isize);
            *fresh0 += 1 as i32 as f32
        }
        j = j.wrapping_add(1)
    };
}
/* * compute the histogram ignoring null elements */
#[no_mangle]
pub unsafe extern "C" fn aubio_hist_do_notnull(s: *mut aubio_hist_t,
                                               input: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut tmp: sint_t = 0 as i32;
    aubio_scale_do((*s).scaler, input);
    /* reset data */
    fvec_zeros((*s).hist);
    /* run accum */
    j = 0 as i32 as uint_t;
    while j < (*input).length {
        if *(*input).data.offset(j as isize) !=
               0 as i32 as f32 {
            tmp = floorf(*(*input).data.offset(j as isize)) as sint_t;
            if tmp >= 0 as i32 && tmp < (*s).nelems as sint_t {
                let ref mut fresh1 = *(*(*s).hist).data.offset(tmp as isize);
                *fresh1 += 1 as i32 as f32
            }
        }
        j = j.wrapping_add(1)
    };
}
/* * compute dynamic histogram for non-null elements */
#[no_mangle]
pub unsafe extern "C" fn aubio_hist_dyn_notnull(s: *mut aubio_hist_t,
                                                input: *mut fvec_t) {
    let mut i: uint_t = 0;
    let mut tmp: sint_t = 0 as i32;
    let ilow: smpl_t = fvec_min(input);
    let ihig: smpl_t = fvec_max(input);
    let step: smpl_t = (ihig - ilow) / (*s).nelems as smpl_t;
    /* readapt */
    aubio_scale_set_limits((*s).scaler, ilow, ihig,
                           0 as i32 as smpl_t, (*s).nelems as smpl_t);
    /* recalculate centers */
    *(*(*s).cent).data.offset(0 as i32 as isize) =
        ilow + 0.5f32 * step;
    i = 1 as i32 as uint_t;
    while i < (*s).nelems {
        *(*(*s).cent).data.offset(i as isize) =
            *(*(*s).cent).data.offset(0 as i32 as isize) +
                i as f32 * step;
        i = i.wrapping_add(1)
    }
    /* scale */
    aubio_scale_do((*s).scaler, input);
    /* reset data */
    fvec_zeros((*s).hist);
    /* run accum */
    i = 0 as i32 as uint_t;
    while i < (*input).length {
        if *(*input).data.offset(i as isize) !=
               0 as i32 as f32 {
            tmp = floorf(*(*input).data.offset(i as isize)) as sint_t;
            if tmp >= 0 as i32 && tmp < (*s).nelems as sint_t {
                let ref mut fresh2 = *(*(*s).hist).data.offset(tmp as isize);
                *fresh2 += 1 as i32 as f32
            }
        }
        i = i.wrapping_add(1)
    };
}
/* * weight the histogram */
#[no_mangle]
pub unsafe extern "C" fn aubio_hist_weight(s: *mut aubio_hist_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*s).nelems {
        let ref mut fresh3 = *(*(*s).hist).data.offset(j as isize);
        *fresh3 *= *(*(*s).cent).data.offset(j as isize);
        j = j.wrapping_add(1)
    };
}
/* * compute the mean of the histogram */
#[no_mangle]
pub unsafe extern "C" fn aubio_hist_mean(s: *const aubio_hist_t)
 -> smpl_t {
    let mut j: uint_t = 0;
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*s).nelems {
        tmp += *(*(*s).hist).data.offset(j as isize);
        j = j.wrapping_add(1)
    }
    return tmp / (*s).nelems as smpl_t;
}
