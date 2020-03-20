extern "C" {
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
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
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_scale_t {
    pub ilow: smpl_t,
    pub ihig: smpl_t,
    pub olow: smpl_t,
    pub ohig: smpl_t,
    pub scaler: smpl_t,
    pub irange: smpl_t,
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
pub type aubio_scale_t = _aubio_scale_t;
/* * create a scale object

  \param flow lower value of output function
  \param fhig higher value of output function
  \param ilow lower value of input function
  \param ihig higher value of output function

*/
/* not implemented yet : type in/out data
     bool inint;
     bool outint;
     */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_scale(mut ilow: smpl_t, mut ihig: smpl_t,
                                         mut olow: smpl_t, mut ohig: smpl_t)
 -> *mut aubio_scale_t {
    let mut s: *mut aubio_scale_t =
        calloc(::std::mem::size_of::<aubio_scale_t>() as u64,
               1 as i32 as u64) as *mut aubio_scale_t;
    aubio_scale_set_limits(s, ilow, ihig, olow, ohig);
    return s;
}
/* * delete a scale object

  \param s scale object as returned by new_aubio_scale

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_scale(mut s: *mut aubio_scale_t) {
    free(s as *mut core::ffi::c_void);
}
/* * modify scale parameters after object creation

  \param s scale object as returned by new_aubio_scale
  \param olow lower value of output function
  \param ohig higher value of output function
  \param ilow lower value of input function
  \param ihig higher value of output function

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_scale_set_limits(mut s: *mut aubio_scale_t,
                                                mut ilow: smpl_t,
                                                mut ihig: smpl_t,
                                                mut olow: smpl_t,
                                                mut ohig: smpl_t) -> uint_t {
    let mut inputrange: smpl_t = ihig - ilow;
    let mut outputrange: smpl_t = ohig - olow;
    (*s).ilow = ilow;
    (*s).ihig = ihig;
    (*s).olow = olow;
    (*s).ohig = ohig;
    if inputrange == 0 as i32 as f32 {
        (*s).scaler = 0.0f64 as smpl_t
    } else {
        (*s).scaler = outputrange / inputrange;
        if inputrange < 0 as i32 as f32 {
            inputrange = inputrange * -1.0f32
        }
    }
    return AUBIO_OK as i32 as uint_t;
}
/* * scale input vector

  \param s scale object as returned by new_aubio_scale
  \param input vector to scale

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_scale_do(mut s: *mut aubio_scale_t,
                                        mut input: *mut fvec_t) {
    let mut j: uint_t = 0;
    j = 0 as i32 as uint_t;
    while j < (*input).length {
        let ref mut fresh0 = *(*input).data.offset(j as isize);
        *fresh0 -= (*s).ilow;
        let ref mut fresh1 = *(*input).data.offset(j as isize);
        *fresh1 *= (*s).scaler;
        let ref mut fresh2 = *(*input).data.offset(j as isize);
        *fresh2 += (*s).olow;
        j = j.wrapping_add(1)
    };
}
