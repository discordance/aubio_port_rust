use ::libc;
extern "C" {
    pub type _aubio_fft_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
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
    /* * revert order of vector elements

  \param s vector to revert

*/
    #[no_mangle]
    fn fvec_rev(s: *mut fvec_t);
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
    /* * make a copy of a vector, applying weights to each element

  \param in input vector
  \param weight weights vector
  \param out output vector

*/
    #[no_mangle]
    fn fvec_weighted_copy(in_0: *const fvec_t, weight: *const fvec_t,
                          out: *mut fvec_t);
    /* * find the index of the min of a vector

  \param s vector to get the index from

  \return the index of the minimum element of v

*/
    #[no_mangle]
    fn fvec_min_elem(s: *mut fvec_t) -> uint_t;
    /* * compute the sum of all elements of a vector

  \param v vector to compute the sum of

  \return the sum of v

*/
    #[no_mangle]
    fn fvec_sum(v: *mut fvec_t) -> smpl_t;
    /* * add a constant to each elements of a vector

  \param v vector to add constant to
  \param c constant to add to v

*/
    #[no_mangle]
    fn fvec_add(v: *mut fvec_t, c: smpl_t);
    /* * finds exact peak index by quadratic interpolation

  See [Quadratic Interpolation of Spectral
  Peaks](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Peak_Interpolation.html),
  by Julius O. Smith III

  \f$ p_{frac} = \frac{1}{2} \frac {x[p-1] - x[p+1]} {x[p-1] - 2 x[p] + x[p+1]} \in [ -.5, .5] \f$

  \param x vector to get the interpolated peak position from
  \param p index of the peak in vector `x`
  \return \f$ p + p_{frac} \f$ exact peak position of interpolated maximum or minimum

*/
    #[no_mangle]
    fn fvec_quadratic_peak_pos(x: *const fvec_t, p: uint_t) -> smpl_t;
    /* * create new FFT computation object

  \param size length of the FFT

*/
    #[no_mangle]
    fn new_aubio_fft(size: uint_t) -> *mut aubio_fft_t;
    /* * delete FFT object

  \param s fft object as returned by new_aubio_fft

*/
    #[no_mangle]
    fn del_aubio_fft(s: *mut aubio_fft_t);
    /* * compute forward FFT

  \param s fft object as returned by new_aubio_fft
  \param input real input signal
  \param compspec complex output fft real/imag

*/
    #[no_mangle]
    fn aubio_fft_do_complex(s: *mut aubio_fft_t, input: *const fvec_t,
                            compspec: *mut fvec_t);
    /* * compute backward (inverse) FFT from real/imag

  \param s fft object as returned by new_aubio_fft
  \param compspec real/imag input fft array
  \param output real output array

*/
    #[no_mangle]
    fn aubio_fft_rdo_complex(s: *mut aubio_fft_t, compspec: *const fvec_t,
                             output: *mut fvec_t);
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
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
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

  Fast Fourier Transform

  Depending on how aubio was compiled, FFT are computed using one of:
    - [Ooura](http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)
    - [FFTW3](http://www.fftw.org)
    - [vDSP](https://developer.apple.com/library/mac/#documentation/Accelerate/Reference/vDSPRef/Reference/reference.html)

  \example spectral/test-fft.c

*/
/* * FFT object

  This object computes forward and backward FFTs.

*/
pub type aubio_fft_t = _aubio_fft_t;
/*
  Copyright (C) 2003-2017 Paul Brossier <piem@aubio.org>

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
/* This algorithm was developed by A. de Cheveigné and H. Kawahara and
 * published in:
 *
 * de Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
 * estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.
 *
 * see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html
 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchyinfast_t {
    pub yin: *mut fvec_t,
    pub tol: smpl_t,
    pub peak_pos: uint_t,
    pub tmpdata: *mut fvec_t,
    pub sqdiff: *mut fvec_t,
    pub kernel: *mut fvec_t,
    pub samples_fft: *mut fvec_t,
    pub kernel_fft: *mut fvec_t,
    pub fft: *mut aubio_fft_t,
}
/*
  Copyright (C) 2003-2017 Paul Brossier <piem@aubio.org>

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

  Pitch detection using YIN algorithm (fast implementation)

  This algorithm was developed by A. de Cheveigne and H. Kawahara and
  published in:

  De Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
  estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.

  This implementation compute the autocorrelation function using time domain
  convolution computed in the spectral domain.

  see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html
      http://recherche.ircam.fr/equipes/pcm/cheveign/ps/2002_JASA_YIN_proof.pdf

*/
/* * pitch detection object */
pub type aubio_pitchyinfast_t = _aubio_pitchyinfast_t;
/* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchyinfast(mut bufsize: uint_t)
 -> *mut aubio_pitchyinfast_t {
    let mut o: *mut aubio_pitchyinfast_t =
        calloc(::std::mem::size_of::<aubio_pitchyinfast_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as
            *mut aubio_pitchyinfast_t;
    (*o).yin =
        new_fvec(bufsize.wrapping_div(2 as libc::c_int as libc::c_uint));
    (*o).tmpdata = new_fvec(bufsize);
    (*o).sqdiff =
        new_fvec(bufsize.wrapping_div(2 as libc::c_int as libc::c_uint));
    (*o).kernel = new_fvec(bufsize);
    (*o).samples_fft = new_fvec(bufsize);
    (*o).kernel_fft = new_fvec(bufsize);
    (*o).fft = new_aubio_fft(bufsize);
    if (*o).yin.is_null() || (*o).tmpdata.is_null() || (*o).tmpdata.is_null()
           || (*o).sqdiff.is_null() || (*o).kernel.is_null() ||
           (*o).samples_fft.is_null() || (*o).kernel.is_null() ||
           (*o).fft.is_null() {
        del_aubio_pitchyinfast(o);
        return 0 as *mut aubio_pitchyinfast_t
    }
    (*o).tol = 0.15f64 as smpl_t;
    (*o).peak_pos = 0 as libc::c_int as uint_t;
    return o;
}
/* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchyin()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchyinfast(mut o:
                                                    *mut aubio_pitchyinfast_t) {
    if !(*o).yin.is_null() { del_fvec((*o).yin); }
    if !(*o).tmpdata.is_null() { del_fvec((*o).tmpdata); }
    if !(*o).sqdiff.is_null() { del_fvec((*o).sqdiff); }
    if !(*o).kernel.is_null() { del_fvec((*o).kernel); }
    if !(*o).samples_fft.is_null() { del_fvec((*o).samples_fft); }
    if !(*o).kernel_fft.is_null() { del_fvec((*o).kernel_fft); }
    if !(*o).fft.is_null() { del_aubio_fft((*o).fft); }
    free(o as *mut libc::c_void);
}
/* * execute pitch detection an input buffer

  \param o pitch detection object as returned by new_aubio_pitchyin()
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
/* all the above in one */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfast_do(mut o:
                                                   *mut aubio_pitchyinfast_t,
                                               mut input: *const fvec_t,
                                               mut out: *mut fvec_t) {
    let tol: smpl_t = (*o).tol; // B / 2
    let mut yin: *mut fvec_t = (*o).yin;
    let length: uint_t = (*yin).length;
    let mut B: uint_t = (*(*o).tmpdata).length;
    let mut W: uint_t = (*(*o).yin).length;
    let mut tmp_slice: fvec_t = fvec_t{length: 0, data: 0 as *mut smpl_t,};
    let mut kernel_ptr: fvec_t = fvec_t{length: 0, data: 0 as *mut smpl_t,};
    let mut tau: uint_t = 0;
    let mut period: sint_t = 0;
    let mut tmp2: smpl_t = 0.0f64 as smpl_t;
    // compute r_t(0) + r_t+tau(0)
    let mut squares: *mut fvec_t = (*o).tmpdata;
    fvec_weighted_copy(input, input, squares);
    tmp_slice.data = (*squares).data;
    tmp_slice.length = W;
    *(*(*o).sqdiff).data.offset(0 as libc::c_int as isize) =
        fvec_sum(&mut tmp_slice);
    tau = 1 as libc::c_int as uint_t;
    while tau < W {
        *(*(*o).sqdiff).data.offset(tau as isize) =
            *(*(*o).sqdiff).data.offset(tau.wrapping_sub(1 as libc::c_int as
                                                             libc::c_uint) as
                                            isize);
        let ref mut fresh0 = *(*(*o).sqdiff).data.offset(tau as isize);
        *fresh0 -=
            *(*squares).data.offset(tau.wrapping_sub(1 as libc::c_int as
                                                         libc::c_uint) as
                                        isize);
        let ref mut fresh1 = *(*(*o).sqdiff).data.offset(tau as isize);
        *fresh1 +=
            *(*squares).data.offset(W.wrapping_add(tau).wrapping_sub(1 as
                                                                         libc::c_int
                                                                         as
                                                                         libc::c_uint)
                                        as isize);
        tau = tau.wrapping_add(1)
    }
    fvec_add((*o).sqdiff,
             *(*(*o).sqdiff).data.offset(0 as libc::c_int as isize));
    // compute r_t(tau) = -2.*ifft(fft(samples)*fft(samples[W-1::-1]))
    let mut compmul: *mut fvec_t = (*o).tmpdata;
    let mut rt_of_tau: *mut fvec_t = (*o).samples_fft;
    aubio_fft_do_complex((*o).fft, input, (*o).samples_fft);
    // build kernel, take a copy of first half of samples
    tmp_slice.data = (*input).data;
    tmp_slice.length = W;
    kernel_ptr.data = (*(*o).kernel).data.offset(1 as libc::c_int as isize);
    kernel_ptr.length = W;
    fvec_copy(&mut tmp_slice, &mut kernel_ptr);
    // reverse them
    fvec_rev(&mut kernel_ptr);
    // compute fft(kernel)
    aubio_fft_do_complex((*o).fft, (*o).kernel, (*o).kernel_fft);
    // compute complex product
    *(*compmul).data.offset(0 as libc::c_int as isize) =
        *(*(*o).kernel_fft).data.offset(0 as libc::c_int as isize) *
            *(*(*o).samples_fft).data.offset(0 as libc::c_int as isize);
    tau = 1 as libc::c_int as uint_t;
    while tau < W {
        *(*compmul).data.offset(tau as isize) =
            *(*(*o).kernel_fft).data.offset(tau as isize) *
                *(*(*o).samples_fft).data.offset(tau as isize);
        let ref mut fresh2 = *(*compmul).data.offset(tau as isize);
        *fresh2 -=
            *(*(*o).kernel_fft).data.offset(B.wrapping_sub(tau) as isize) *
                *(*(*o).samples_fft).data.offset(B.wrapping_sub(tau) as
                                                     isize);
        tau = tau.wrapping_add(1)
    }
    *(*compmul).data.offset(W as isize) =
        *(*(*o).kernel_fft).data.offset(W as isize) *
            *(*(*o).samples_fft).data.offset(W as isize);
    tau = 1 as libc::c_int as uint_t;
    while tau < W {
        *(*compmul).data.offset(B.wrapping_sub(tau) as isize) =
            *(*(*o).kernel_fft).data.offset(B.wrapping_sub(tau) as isize) *
                *(*(*o).samples_fft).data.offset(tau as isize);
        let ref mut fresh3 =
            *(*compmul).data.offset(B.wrapping_sub(tau) as isize);
        *fresh3 +=
            *(*(*o).kernel_fft).data.offset(tau as isize) *
                *(*(*o).samples_fft).data.offset(B.wrapping_sub(tau) as
                                                     isize);
        tau = tau.wrapping_add(1)
    }
    // compute inverse fft
    aubio_fft_rdo_complex((*o).fft, compmul, rt_of_tau);
    // compute square difference r_t(tau) = sqdiff - 2 * r_t_tau[W-1:-1]
    tau = 0 as libc::c_int as uint_t;
    while tau < W {
        *(*yin).data.offset(tau as isize) =
            (*(*(*o).sqdiff).data.offset(tau as isize) as libc::c_double -
                 2.0f64 *
                     *(*rt_of_tau).data.offset(tau.wrapping_add(W) as isize)
                         as libc::c_double) as smpl_t;
        tau = tau.wrapping_add(1)
    }
    // now build yin and look for first minimum
    fvec_zeros(out);
    *(*yin).data.offset(0 as libc::c_int as isize) = 1.0f64 as smpl_t;
    tau = 1 as libc::c_int as uint_t;
    while tau < length {
        tmp2 += *(*yin).data.offset(tau as isize);
        if tmp2 != 0 as libc::c_int as libc::c_float {
            let ref mut fresh4 = *(*yin).data.offset(tau as isize);
            *fresh4 *= tau as libc::c_float / tmp2
        } else { *(*yin).data.offset(tau as isize) = 1.0f64 as smpl_t }
        period = tau.wrapping_sub(3 as libc::c_int as libc::c_uint) as sint_t;
        if tau > 4 as libc::c_int as libc::c_uint &&
               *(*yin).data.offset(period as isize) < tol &&
               *(*yin).data.offset(period as isize) <
                   *(*yin).data.offset((period + 1 as libc::c_int) as isize) {
            (*o).peak_pos = period as uint_t;
            *(*out).data.offset(0 as libc::c_int as isize) =
                fvec_quadratic_peak_pos(yin, (*o).peak_pos);
            return
        }
        tau = tau.wrapping_add(1)
    }
    // use global minimum 
    (*o).peak_pos = fvec_min_elem(yin);
    *(*out).data.offset(0 as libc::c_int as isize) =
        fvec_quadratic_peak_pos(yin, (*o).peak_pos);
}
/* * get current confidence of YIN algorithm

  \param o YIN pitch detection object
  \return confidence parameter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfast_get_confidence(mut o:
                                                               *mut aubio_pitchyinfast_t)
 -> smpl_t {
    return (1.0f64 -
                *(*(*o).yin).data.offset((*o).peak_pos as isize) as
                    libc::c_double) as smpl_t;
}
/* * set tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \param tol tolerance parameter for minima selection [default 0.15]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfast_set_tolerance(mut o:
                                                              *mut aubio_pitchyinfast_t,
                                                          mut tol: smpl_t)
 -> uint_t {
    (*o).tol = tol;
    return 0 as libc::c_int as uint_t;
}
/* * get tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \return tolerance parameter for minima selection [default 0.15]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfast_get_tolerance(mut o:
                                                              *mut aubio_pitchyinfast_t)
 -> smpl_t {
    return (*o).tol;
}
