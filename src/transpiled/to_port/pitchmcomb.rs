use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn fabsf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn powf(_: libc::c_float, _: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
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
    /* * compute the principal argument

  This function maps the input phase to its corresponding value wrapped in the
range \f$ [-\pi, \pi] \f$.

  \param phase unwrapped phase to map to the unit circle

  \return equivalent phase wrapped to the unit circle

*/
    #[no_mangle]
    fn aubio_unwrap2pi(phase: smpl_t) -> smpl_t;
    /* *  alpha normalisation

  This function divides all elements of a vector by the p-norm as computed by
fvec_alpha_norm().

  \param v vector to compute norm from
  \param p order of the computed norm

*/
    #[no_mangle]
    fn fvec_alpha_normalise(v: *mut fvec_t, p: smpl_t);
    /* * add a constant to each elements of a vector

  \param v vector to add constant to
  \param c constant to add to v

*/
    #[no_mangle]
    fn fvec_add(v: *mut fvec_t, c: smpl_t);
    /* * remove the minimum value of the vector to each elements

  \param v vector to remove minimum from

*/
    #[no_mangle]
    fn fvec_min_removal(v: *mut fvec_t);
    /* * apply adaptive threshold to a vector

  For each points at position p of an input vector, this function remove the
moving median threshold computed at p.

  \param v input vector
  \param tmp temporary vector of length post+1+pre
  \param post length of causal part to take before pos
  \param pre length of anti-causal part to take after pos

*/
    #[no_mangle]
    fn fvec_adapt_thres(v: *mut fvec_t, tmp: *mut fvec_t, post: uint_t,
                        pre: uint_t);
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
    #[no_mangle]
    fn fvec_peakpick(v: *const fvec_t, p: uint_t) -> uint_t;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchmcomb_t {
    pub threshold: smpl_t,
    pub alpha: smpl_t,
    pub cutoff: smpl_t,
    pub tol: smpl_t,
    pub win_post: uint_t,
    pub win_pre: uint_t,
    pub ncand: uint_t,
    pub npartials: uint_t,
    pub count: uint_t,
    pub goodcandidate: uint_t,
    pub spec_partition: uint_t,
    pub peaks: *mut aubio_spectralpeak_t,
    pub candidates: *mut *mut aubio_spectralcandidate_t,
    pub newmag: *mut fvec_t,
    pub scratch: *mut fvec_t,
    pub scratch2: *mut fvec_t,
    pub theta: *mut fvec_t,
    pub phasediff: smpl_t,
    pub phasefreq: smpl_t,
}
pub type aubio_spectralcandidate_t = _aubio_spectralcandidate_t;
/* *< peak magnitude */
/* * spectral candidates array object */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_spectralcandidate_t {
    pub ebin: smpl_t,
    pub ecomb: *mut smpl_t,
    pub ene: smpl_t,
    pub len: smpl_t,
}
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
pub type aubio_spectralpeak_t = _aubio_spectralpeak_t;
/* * threshfn: name or handle of fn for computing adaptive threshold [median] */
  /* * aubio_thresholdfn_t thresholdfn; */
  /* * picker: name or handle of fn for picking event times [quadpick] */
  /* * aubio_pickerfn_t pickerfn; */
/* * spectral peak object */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_spectralpeak_t {
    pub bin: uint_t,
    pub ebin: smpl_t,
    pub mag: smpl_t,
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

  Pitch detection using multiple-comb filter

  This fundamental frequency estimation algorithm implements spectral
  flattening, multi-comb filtering and peak histogramming.

  This method was designed by Juan P. Bello and described in:

  Juan-Pablo Bello. ``Towards the Automated Analysis of Simple Polyphonic
  Music''.  PhD thesis, Centre for Digital Music, Queen Mary University of
  London, London, UK, 2003.

  \example pitch/test-pitchmcomb.c

*/
/* * pitch detection object */
pub type aubio_pitchmcomb_t = _aubio_pitchmcomb_t;
/* * execute pitch detection on an input spectral frame

  \param p pitch detection object as returned by new_aubio_pitchmcomb
  \param in_fftgrain input signal spectrum as computed by aubio_pvoc_do
  \param out_cands pitch candidate frequenciess, in bins

*/
/* *< length */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchmcomb_do(mut p: *mut aubio_pitchmcomb_t,
                                             mut fftgrain: *const cvec_t,
                                             mut output: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut instfreq: smpl_t = 0.;
    let mut newmag: *mut fvec_t = (*p).newmag;
    //smpl_t hfc; //fe=instfreq(theta1,theta,ops); //theta1=theta;
  /* copy incoming grain to newmag */
    j = 0 as libc::c_int as uint_t;
    while j < (*newmag).length {
        *(*newmag).data.offset(j as isize) =
            *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
    /* detect only if local energy > 10. */
  //if (aubio_level_lin (newmag) * newmag->length > 10.) {
  //hfc = fvec_local_hfc(newmag); //not used
    aubio_pitchmcomb_spectral_pp(p, newmag);
    aubio_pitchmcomb_combdet(p, newmag);
    //aubio_pitchmcomb_sort_cand_freq(p->candidates,p->ncand);
  //return p->candidates[p->goodcandidate]->ebin;
    j =
        floorf(((**(*p).candidates.offset((*p).goodcandidate as isize)).ebin
                    as libc::c_double + 0.5f64) as libc::c_float) as uint_t;
    instfreq =
        aubio_unwrap2pi(*(*fftgrain).phas.offset(j as isize) -
                            *(*(*p).theta).data.offset(j as isize) -
                            j as libc::c_float * (*p).phasediff);
    instfreq *= (*p).phasefreq;
    /* store phase for next run */
    j = 0 as libc::c_int as uint_t;
    while j < (*(*p).theta).length {
        *(*(*p).theta).data.offset(j as isize) =
            *(*fftgrain).phas.offset(j as isize);
        j = j.wrapping_add(1)
    }
    //return p->candidates[p->goodcandidate]->ebin;
    *(*output).data.offset(0 as libc::c_int as isize) =
        floorf(((**(*p).candidates.offset((*p).goodcandidate as isize)).ebin
                    as libc::c_double + 0.5f64) as libc::c_float) + instfreq;
    /*} else {
     return -1.;
     } */
}
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchmcomb_spectral_pp(mut p:
                                                          *mut aubio_pitchmcomb_t,
                                                      mut newmag:
                                                          *const fvec_t) {
    let mut mag: *mut fvec_t = (*p).scratch;
    let mut tmp: *mut fvec_t = (*p).scratch2;
    let mut j: uint_t = 0;
    let mut length: uint_t = (*mag).length;
    /* copy newmag to mag (scracth) */
    j = 0 as libc::c_int as uint_t; /* min removal          */
    while j < length {
        *(*mag).data.offset(j as isize) =
            *(*newmag).data.offset(j as isize); /* alpha normalisation  */
        j = j.wrapping_add(1)
    }
    fvec_min_removal(mag);
    fvec_alpha_normalise(mag, (*p).alpha);
    /* skipped */
    /* low pass filtering   */
    /* * \bug fvec_moving_thres may write out of bounds */
    fvec_adapt_thres(mag, tmp, (*p).win_post,
                     (*p).win_pre); /* adaptative threshold */
    fvec_add(mag, -(*p).threshold); /* fixed threshold      */
    let mut peaks: *mut aubio_spectralpeak_t = (*p).peaks;
    let mut count: uint_t = 0;
    /*  return bin and ebin */
    count = aubio_pitchmcomb_quadpick(peaks, mag);
    j = 0 as libc::c_int as uint_t;
    while j < count {
        (*peaks.offset(j as isize)).mag =
            *(*newmag).data.offset((*peaks.offset(j as isize)).bin as isize);
        j = j.wrapping_add(1)
    }
    /* reset non peaks */
    j = count;
    while j < length {
        (*peaks.offset(j as isize)).mag = 0.0f64 as smpl_t;
        j = j.wrapping_add(1)
    }
    (*p).peaks = peaks;
    (*p).count = count;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchmcomb_combdet(mut p:
                                                      *mut aubio_pitchmcomb_t,
                                                  mut newmag: *const fvec_t) {
    let mut peaks: *mut aubio_spectralpeak_t = (*p).peaks;
    let mut candidate: *mut *mut aubio_spectralcandidate_t = (*p).candidates;
    /* parms */
    let mut N: uint_t =
        (*p).npartials; /* maximum number of partials to be considered 10 */
    let mut M: uint_t =
        (*p).ncand; /* maximum number of combs to be considered 5 */
    let mut length: uint_t = (*newmag).length;
    let mut count: uint_t = (*p).count;
    let mut k: uint_t = 0;
    let mut l: uint_t = 0;
    let mut d: uint_t = 0;
    let mut curlen: uint_t = 0 as libc::c_int as uint_t;
    let mut delta2: smpl_t = 0.;
    let mut xx: smpl_t = 0.;
    let mut position: uint_t = 0 as libc::c_int as uint_t;
    let mut root_peak: uint_t = 0 as libc::c_int as uint_t;
    let mut tmpl: uint_t = 0 as libc::c_int as uint_t;
    let mut tmpene: smpl_t = 0.0f64 as smpl_t;
    /* get the biggest peak in the spectrum */
    root_peak = aubio_pitchmcomb_get_root_peak(peaks, count);
    /* not enough partials in highest notes, could be forced */
  //if (peaks[root_peak].ebin >= aubio_miditofreq(85.)/p->tau) N=2;
  //if (peaks[root_peak].ebin >= aubio_miditofreq(90.)/p->tau) N=1;
  /* now calculate the energy of each of the 5 combs */
    l = 0 as libc::c_int as uint_t; /* reset ene and len sums */
    while l < M {
        let mut scaler: smpl_t =
            (1.0f64 / (l as libc::c_double + 1.0f64)) as smpl_t;
        (**candidate.offset(l as isize)).ene = 0.0f64 as smpl_t;
        (**candidate.offset(l as isize)).len = 0.0f64 as smpl_t;
        (**candidate.offset(l as isize)).ebin =
            scaler * (*peaks.offset(root_peak as isize)).ebin;
        /* if less than N peaks available, curlen < N */
        if (**candidate.offset(l as isize)).ebin as libc::c_double != 0.0f64 {
            curlen =
                floorf(length as libc::c_float /
                           (**candidate.offset(l as isize)).ebin) as uint_t
        }
        curlen = if N < curlen { N } else { curlen };
        /* fill candidate[l]->ecomb[k] with (k+1)*candidate[l]->ebin */
        k = 0 as libc::c_int as uint_t;
        while k < curlen {
            *(**candidate.offset(l as isize)).ecomb.offset(k as isize) =
                ((**candidate.offset(l as isize)).ebin as libc::c_double *
                     (k as libc::c_double + 1.0f64)) as smpl_t;
            k = k.wrapping_add(1)
        }
        k = curlen;
        while k < length {
            *(**candidate.offset(l as isize)).ecomb.offset(k as isize) =
                0.0f64 as smpl_t;
            k = k.wrapping_add(1)
        }
        /* for each in candidate[l]->ecomb[k] */
        k = 0 as libc::c_int as uint_t;
        while k < curlen {
            xx = 100000.0f64 as smpl_t;
            /* * get the candidate->ecomb the closer to peaks.ebin
       * (to cope with the inharmonicity)*/
            d = 0 as libc::c_int as uint_t;
            while d < count {
                delta2 =
                    fabsf(*(**candidate.offset(l as
                                                   isize)).ecomb.offset(k as
                                                                            isize)
                              - (*peaks.offset(d as isize)).ebin);
                if delta2 <= xx { position = d; xx = delta2 }
                d = d.wrapping_add(1)
            }
            /* for a Q factor of 17, maintaining "constant Q filtering",
       * and sum energy and length over non null combs */
            if (17.0f64 * xx as libc::c_double) <
                   *(**candidate.offset(l as isize)).ecomb.offset(k as isize)
                       as libc::c_double {
                *(**candidate.offset(l as isize)).ecomb.offset(k as isize) =
                    (*peaks.offset(position as isize)).ebin;
                let ref mut fresh0 = (**candidate.offset(l as isize)).ene;
                *fresh0 +=
                    powf(*(*newmag).data.offset(floorf((*(**candidate.offset(l
                                                                                 as
                                                                                 isize)).ecomb.offset(k
                                                                                                          as
                                                                                                          isize)
                                                            as libc::c_double
                                                            + 0.5f64) as
                                                           libc::c_float) as
                                                    uint_t as isize),
                         0.25f64 as libc::c_float);
                let ref mut fresh1 = (**candidate.offset(l as isize)).len;
                *fresh1 =
                    (*fresh1 as libc::c_double +
                         1.0f64 / curlen as libc::c_double) as smpl_t
            } else {
                *(**candidate.offset(l as isize)).ecomb.offset(k as isize) =
                    0.0f64 as smpl_t
            }
            k = k.wrapping_add(1)
        }
        /* punishment */
    /*if (candidate[l]->len<0.6)
       candidate[l]->ene=0.; */
    /* remember best candidate energy (in polyphonic, could check for
     * tmpene*1.1 < candidate->ene to reduce jumps towards low frequencies) */
        if tmpene < (**candidate.offset(l as isize)).ene {
            tmpl = l;
            tmpene = (**candidate.offset(l as isize)).ene
        }
        l = l.wrapping_add(1)
    }
    //p->candidates=candidate;
  //p->peaks=peaks;
    (*p).goodcandidate = tmpl;
}
/* * T=quadpick(X): return indices of elements of X which are peaks and positive
 * exact peak positions are retrieved by quadratic interpolation
 *
 * \bug peak-picking too picky, sometimes counts too many peaks ?
 */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchmcomb_quadpick(mut spectral_peaks:
                                                       *mut aubio_spectralpeak_t,
                                                   mut X: *const fvec_t)
 -> uint_t {
    let mut j: uint_t = 0;
    let mut ispeak: uint_t = 0;
    let mut count: uint_t = 0 as libc::c_int as uint_t;
    j = 1 as libc::c_int as uint_t;
    while j < (*X).length.wrapping_sub(1 as libc::c_int as libc::c_uint) {
        ispeak = fvec_peakpick(X, j);
        if ispeak != 0 {
            count =
                (count as libc::c_uint).wrapping_add(ispeak) as uint_t as
                    uint_t;
            (*spectral_peaks.offset(count.wrapping_sub(1 as libc::c_int as
                                                           libc::c_uint) as
                                        isize)).bin = j;
            (*spectral_peaks.offset(count.wrapping_sub(1 as libc::c_int as
                                                           libc::c_uint) as
                                        isize)).ebin =
                fvec_quadratic_peak_pos(X, j)
        }
        j = j.wrapping_add(1)
    }
    return count;
}
/* get predominant partial */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchmcomb_get_root_peak(mut peaks:
                                                            *mut aubio_spectralpeak_t,
                                                        mut length: uint_t)
 -> uint_t {
    let mut i: uint_t = 0;
    let mut pos: uint_t = 0 as libc::c_int as uint_t;
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    i = 0 as libc::c_int as uint_t;
    while i < length {
        if tmp <= (*peaks.offset(i as isize)).mag {
            pos = i;
            tmp = (*peaks.offset(i as isize)).mag
        }
        i = i.wrapping_add(1)
    }
    return pos;
}
/* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse
  \param hop_size step size between two consecutive analysis instant

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchmcomb(mut bufsize: uint_t,
                                              mut hopsize: uint_t)
 -> *mut aubio_pitchmcomb_t {
    let mut p: *mut aubio_pitchmcomb_t =
        calloc(::std::mem::size_of::<aubio_pitchmcomb_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_pitchmcomb_t;
    /* bug: should check if size / 8 > post+pre+1 */
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    let mut spec_size: uint_t = 0;
    (*p).spec_partition = 2 as libc::c_int as uint_t;
    (*p).ncand = 5 as libc::c_int as uint_t;
    (*p).npartials = 5 as libc::c_int as uint_t;
    (*p).cutoff = 1.0f64 as smpl_t;
    (*p).threshold = 0.01f64 as smpl_t;
    (*p).win_post = 8 as libc::c_int as uint_t;
    (*p).win_pre = 7 as libc::c_int as uint_t;
    // p->tau              = samplerate/bufsize;
    (*p).alpha = 9.0f64 as smpl_t;
    (*p).goodcandidate = 0 as libc::c_int as uint_t;
    (*p).phasefreq =
        (bufsize.wrapping_div(hopsize) as libc::c_double /
             (3.14159265358979323846264338327950288f64 * 2.0f64)) as smpl_t;
    (*p).phasediff =
        (3.14159265358979323846264338327950288f64 * 2.0f64 *
             hopsize as libc::c_double / bufsize as libc::c_double) as smpl_t;
    spec_size =
        bufsize.wrapping_div((*p).spec_partition).wrapping_add(1 as
                                                                   libc::c_int
                                                                   as
                                                                   libc::c_uint);
    //p->pickerfn = quadpick;
  //p->biquad = new_biquad(0.1600,0.3200,0.1600, -0.5949, 0.2348);
  /* allocate temp memory */
    (*p).newmag = new_fvec(spec_size);
    /* array for median */
    (*p).scratch = new_fvec(spec_size);
    /* array for phase */
    (*p).theta = new_fvec(spec_size);
    /* array for adaptative threshold */
    (*p).scratch2 =
        new_fvec((*p).win_post.wrapping_add((*p).win_pre).wrapping_add(1 as
                                                                           libc::c_int
                                                                           as
                                                                           libc::c_uint));
    /* array of spectral peaks */
    (*p).peaks =
        calloc((spec_size as
                    libc::c_ulong).wrapping_mul(::std::mem::size_of::<aubio_spectralpeak_t>()
                                                    as libc::c_ulong),
               1 as libc::c_int as libc::c_ulong) as
            *mut aubio_spectralpeak_t;
    i = 0 as libc::c_int as uint_t;
    while i < spec_size {
        (*(*p).peaks.offset(i as isize)).bin = 0.0f64 as uint_t;
        (*(*p).peaks.offset(i as isize)).ebin = 0.0f64 as smpl_t;
        (*(*p).peaks.offset(i as isize)).mag = 0.0f64 as smpl_t;
        i = i.wrapping_add(1)
    }
    /* array of pointers to spectral candidates */
    (*p).candidates =
        calloc(((*p).ncand as
                    libc::c_ulong).wrapping_mul(::std::mem::size_of::<*mut aubio_spectralcandidate_t>()
                                                    as libc::c_ulong),
               1 as libc::c_int as libc::c_ulong) as
            *mut *mut aubio_spectralcandidate_t;
    i = 0 as libc::c_int as uint_t;
    while i < (*p).ncand {
        let ref mut fresh2 = *(*p).candidates.offset(i as isize);
        *fresh2 =
            calloc(::std::mem::size_of::<aubio_spectralcandidate_t>() as
                       libc::c_ulong, 1 as libc::c_int as libc::c_ulong) as
                *mut aubio_spectralcandidate_t;
        let ref mut fresh3 = (**(*p).candidates.offset(i as isize)).ecomb;
        *fresh3 =
            calloc((spec_size as
                        libc::c_ulong).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as libc::c_ulong),
                   1 as libc::c_int as libc::c_ulong) as *mut smpl_t;
        j = 0 as libc::c_int as uint_t;
        while j < spec_size {
            *(**(*p).candidates.offset(i as isize)).ecomb.offset(j as isize) =
                0.0f64 as smpl_t;
            j = j.wrapping_add(1)
        }
        (**(*p).candidates.offset(i as isize)).ene = 0.0f64 as smpl_t;
        (**(*p).candidates.offset(i as isize)).ebin = 0.0f64 as smpl_t;
        (**(*p).candidates.offset(i as isize)).len = 0.0f64 as smpl_t;
        i = i.wrapping_add(1)
    }
    return p;
}
/* * deletion of the pitch detection object

  \param p pitch detection object as returned by new_aubio_pitchfcomb

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchmcomb(mut p:
                                                  *mut aubio_pitchmcomb_t) {
    let mut i: uint_t = 0;
    del_fvec((*p).newmag);
    del_fvec((*p).scratch);
    del_fvec((*p).theta);
    del_fvec((*p).scratch2);
    free((*p).peaks as *mut libc::c_void);
    i = 0 as libc::c_int as uint_t;
    while i < (*p).ncand {
        free((**(*p).candidates.offset(i as isize)).ecomb as
                 *mut libc::c_void);
        free(*(*p).candidates.offset(i as isize) as *mut libc::c_void);
        i = i.wrapping_add(1)
    }
    free((*p).candidates as *mut libc::c_void);
    free(p as *mut libc::c_void);
}
