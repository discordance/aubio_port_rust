extern "C" {
    #[no_mangle]
    fn atanf(_: f32) -> f32;
    #[no_mangle]
    fn cosf(_: f32) -> f32;
    #[no_mangle]
    fn sinf(_: f32) -> f32;
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
// modifications made for aubio:
//  - replace all 'double' with 'smpl_t'
//  - include "aubio_priv.h" (for config.h and types.h)
//  - add missing prototypes
//  - use COS, SIN, and ATAN macros
//  - add cast to (smpl_t) to avoid float conversion warnings
//  - declare initialization as static
//  - prefix public function with aubio_ooura_
/*
Fast Fourier/Cosine/Sine Transform
    dimension   :one
    data length :power of 2
    decimation  :frequency
    radix       :8, 4, 2
    data        :inplace
    table       :use
functions
    cdft: Complex Discrete Fourier Transform
    rdft: Real Discrete Fourier Transform
    ddct: Discrete Cosine Transform
    ddst: Discrete Sine Transform
    dfct: Cosine Transform of RDFT (Real Symmetric DFT)
    dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
function prototypes
    void cdft(int, int, smpl_t *, int *, smpl_t *);
    void rdft(int, int, smpl_t *, int *, smpl_t *);
    void ddct(int, int, smpl_t *, int *, smpl_t *);
    void ddst(int, int, smpl_t *, int *, smpl_t *);
    void dfct(int, smpl_t *, smpl_t *, int *, smpl_t *);
    void dfst(int, smpl_t *, smpl_t *, int *, smpl_t *);


-------- Complex DFT (Discrete Fourier Transform) --------
    [definition]
        <case1>
            X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
        <case2>
            X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
        (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
    [usage]
        <case1>
            ip[0] = 0; // first time only
            cdft(2*n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            cdft(2*n, -1, a, ip, w);
    [parameters]
        2*n            :data length (int)
                        n >= 1, n = power of 2
        a[0...2*n-1]   :input/output data (smpl_t *)
                        input data
                            a[2*j] = Re(x[j]), 
                            a[2*j+1] = Im(x[j]), 0<=j<n
                        output data
                            a[2*k] = Re(X[k]), 
                            a[2*k+1] = Im(X[k]), 0<=k<n
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (smpl_t *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            cdft(2*n, -1, a, ip, w);
        is 
            cdft(2*n, 1, a, ip, w);
            for (j = 0; j <= 2 * n - 1; j++) {
                a[j] *= 1.0 / n;
            }
        .


-------- Real DFT / Inverse of Real DFT --------
    [definition]
        <case1> RDFT
            R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
            I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
        <case2> IRDFT (excluding scale)
            a[k] = (R[0] + R[n/2]*cos(pi*k))/2 + 
                   sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) + 
                   sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            rdft(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            rdft(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (smpl_t *)
                        <case1>
                            output data
                                a[2*k] = R[k], 0<=k<n/2
                                a[2*k+1] = I[k], 0<k<n/2
                                a[1] = R[n/2]
                        <case2>
                            input data
                                a[2*j] = R[j], 0<=j<n/2
                                a[2*j+1] = I[j], 0<j<n/2
                                a[1] = R[n/2]
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (smpl_t *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            rdft(n, 1, a, ip, w);
        is 
            rdft(n, -1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
    [definition]
        <case1> IDCT (excluding scale)
            C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
        <case2> DCT
            C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            ddct(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            ddct(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (smpl_t *)
                        output data
                            a[k] = C[k], 0<=k<n
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/4-1] :cos/sin table (smpl_t *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            ddct(n, -1, a, ip, w);
        is 
            a[0] *= 0.5;
            ddct(n, 1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- DST (Discrete Sine Transform) / Inverse of DST --------
    [definition]
        <case1> IDST (excluding scale)
            S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
        <case2> DST
            S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            ddst(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            ddst(n, -1, a, ip, w);
    [parameters]
        n              :data length (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (smpl_t *)
                        <case1>
                            input data
                                a[j] = A[j], 0<j<n
                                a[0] = A[n]
                            output data
                                a[k] = S[k], 0<=k<n
                        <case2>
                            output data
                                a[k] = S[k], 0<k<n
                                a[0] = S[n]
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/4-1] :cos/sin table (smpl_t *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            ddst(n, -1, a, ip, w);
        is 
            a[0] *= 0.5;
            ddst(n, 1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- Cosine Transform of RDFT (Real Symmetric DFT) --------
    [definition]
        C[k] = sum_j=0^n a[j]*cos(pi*j*k/n), 0<=k<=n
    [usage]
        ip[0] = 0; // first time only
        dfct(n, a, t, ip, w);
    [parameters]
        n              :data length - 1 (int)
                        n >= 2, n = power of 2
        a[0...n]       :input/output data (smpl_t *)
                        output data
                            a[k] = C[k], 0<=k<=n
        t[0...n/2]     :work area (smpl_t *)
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/4)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/4+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/8-1] :cos/sin table (smpl_t *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            a[0] *= 0.5;
            a[n] *= 0.5;
            dfct(n, a, t, ip, w);
        is 
            a[0] *= 0.5;
            a[n] *= 0.5;
            dfct(n, a, t, ip, w);
            for (j = 0; j <= n; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
    [definition]
        S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
    [usage]
        ip[0] = 0; // first time only
        dfst(n, a, t, ip, w);
    [parameters]
        n              :data length + 1 (int)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (smpl_t *)
                        output data
                            a[k] = S[k], 0<k<n
                        (a[0] is used for work area)
        t[0...n/2-1]   :work area (smpl_t *)
        ip[0...*]      :work area for bit reversal (int *)
                        length of ip >= 2+sqrt(n/4)
                        strictly, 
                        length of ip >= 
                            2+(1<<(int)(log(n/4+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/8-1] :cos/sin table (smpl_t *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            dfst(n, a, t, ip, w);
        is 
            dfst(n, a, t, ip, w);
            for (j = 1; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


Appendix :
    The cos/sin table is recalculated when the larger table required.
    w[] and ip[] are compatible with all routines.
*/
#[no_mangle]
pub unsafe extern "C" fn aubio_ooura_cdft(n: i32,
                                          isgn: i32,
                                          a: *mut smpl_t,
                                          ip: *mut i32,
                                          w: *mut smpl_t) {
    if n > *ip.offset(0 as i32 as isize) << 2 as i32 {
        makewt(n >> 2 as i32, ip, w);
    }
    if n > 4 as i32 {
        if isgn >= 0 as i32 {
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftfsub(n, a, w);
        } else {
            bitrv2conj(n, ip.offset(2 as i32 as isize), a);
            cftbsub(n, a, w);
        }
    } else if n == 4 as i32 { cftfsub(n, a, w); };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_ooura_rdft(n: i32,
                                          isgn: i32,
                                          a: *mut smpl_t,
                                          ip: *mut i32,
                                          w: *mut smpl_t) {
    let mut nw: i32 = 0;
    let mut nc: i32 = 0;
    let mut xi: smpl_t = 0.;
    nw = *ip.offset(0 as i32 as isize);
    if n > nw << 2 as i32 {
        nw = n >> 2 as i32;
        makewt(nw, ip, w);
    }
    nc = *ip.offset(1 as i32 as isize);
    if n > nc << 2 as i32 {
        nc = n >> 2 as i32;
        makect(nc, ip, w.offset(nw as isize));
    }
    if isgn >= 0 as i32 {
        if n > 4 as i32 {
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w.offset(nw as isize));
        } else if n == 4 as i32 { cftfsub(n, a, w); }
        xi =
            *a.offset(0 as i32 as isize) -
                *a.offset(1 as i32 as isize);
        let ref mut fresh0 = *a.offset(0 as i32 as isize);
        *fresh0 += *a.offset(1 as i32 as isize);
        *a.offset(1 as i32 as isize) = xi
    } else {
        *a.offset(1 as i32 as isize) =
            0.5f64 as smpl_t *
                (*a.offset(0 as i32 as isize) -
                     *a.offset(1 as i32 as isize));
        let ref mut fresh1 = *a.offset(0 as i32 as isize);
        *fresh1 -= *a.offset(1 as i32 as isize);
        if n > 4 as i32 {
            rftbsub(n, a, nc, w.offset(nw as isize));
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftbsub(n, a, w);
        } else if n == 4 as i32 { cftfsub(n, a, w); }
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_ooura_ddct(n: i32,
                                          isgn: i32,
                                          a: *mut smpl_t,
                                          ip: *mut i32,
                                          w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut nw: i32 = 0;
    let mut nc: i32 = 0;
    let mut xr: smpl_t = 0.;
    nw = *ip.offset(0 as i32 as isize);
    if n > nw << 2 as i32 {
        nw = n >> 2 as i32;
        makewt(nw, ip, w);
    }
    nc = *ip.offset(1 as i32 as isize);
    if n > nc { nc = n; makect(nc, ip, w.offset(nw as isize)); }
    if isgn < 0 as i32 {
        xr = *a.offset((n - 1 as i32) as isize);
        j = n - 2 as i32;
        while j >= 2 as i32 {
            *a.offset((j + 1 as i32) as isize) =
                *a.offset(j as isize) -
                    *a.offset((j - 1 as i32) as isize);
            let ref mut fresh2 = *a.offset(j as isize);
            *fresh2 += *a.offset((j - 1 as i32) as isize);
            j -= 2 as i32
        }
        *a.offset(1 as i32 as isize) =
            *a.offset(0 as i32 as isize) - xr;
        let ref mut fresh3 = *a.offset(0 as i32 as isize);
        *fresh3 += xr;
        if n > 4 as i32 {
            rftbsub(n, a, nc, w.offset(nw as isize));
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftbsub(n, a, w);
        } else if n == 4 as i32 { cftfsub(n, a, w); }
    }
    dctsub(n, a, nc, w.offset(nw as isize));
    if isgn >= 0 as i32 {
        if n > 4 as i32 {
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w.offset(nw as isize));
        } else if n == 4 as i32 { cftfsub(n, a, w); }
        xr =
            *a.offset(0 as i32 as isize) -
                *a.offset(1 as i32 as isize);
        let ref mut fresh4 = *a.offset(0 as i32 as isize);
        *fresh4 += *a.offset(1 as i32 as isize);
        j = 2 as i32;
        while j < n {
            *a.offset((j - 1 as i32) as isize) =
                *a.offset(j as isize) -
                    *a.offset((j + 1 as i32) as isize);
            let ref mut fresh5 = *a.offset(j as isize);
            *fresh5 += *a.offset((j + 1 as i32) as isize);
            j += 2 as i32
        }
        *a.offset((n - 1 as i32) as isize) = xr
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_ooura_ddst(n: i32,
                                          isgn: i32,
                                          a: *mut smpl_t,
                                          ip: *mut i32,
                                          w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut nw: i32 = 0;
    let mut nc: i32 = 0;
    let mut xr: smpl_t = 0.;
    nw = *ip.offset(0 as i32 as isize);
    if n > nw << 2 as i32 {
        nw = n >> 2 as i32;
        makewt(nw, ip, w);
    }
    nc = *ip.offset(1 as i32 as isize);
    if n > nc { nc = n; makect(nc, ip, w.offset(nw as isize)); }
    if isgn < 0 as i32 {
        xr = *a.offset((n - 1 as i32) as isize);
        j = n - 2 as i32;
        while j >= 2 as i32 {
            *a.offset((j + 1 as i32) as isize) =
                -*a.offset(j as isize) -
                    *a.offset((j - 1 as i32) as isize);
            let ref mut fresh6 = *a.offset(j as isize);
            *fresh6 -= *a.offset((j - 1 as i32) as isize);
            j -= 2 as i32
        }
        *a.offset(1 as i32 as isize) =
            *a.offset(0 as i32 as isize) + xr;
        let ref mut fresh7 = *a.offset(0 as i32 as isize);
        *fresh7 -= xr;
        if n > 4 as i32 {
            rftbsub(n, a, nc, w.offset(nw as isize));
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftbsub(n, a, w);
        } else if n == 4 as i32 { cftfsub(n, a, w); }
    }
    dstsub(n, a, nc, w.offset(nw as isize));
    if isgn >= 0 as i32 {
        if n > 4 as i32 {
            bitrv2(n, ip.offset(2 as i32 as isize), a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w.offset(nw as isize));
        } else if n == 4 as i32 { cftfsub(n, a, w); }
        xr =
            *a.offset(0 as i32 as isize) -
                *a.offset(1 as i32 as isize);
        let ref mut fresh8 = *a.offset(0 as i32 as isize);
        *fresh8 += *a.offset(1 as i32 as isize);
        j = 2 as i32;
        while j < n {
            *a.offset((j - 1 as i32) as isize) =
                -*a.offset(j as isize) -
                    *a.offset((j + 1 as i32) as isize);
            let ref mut fresh9 = *a.offset(j as isize);
            *fresh9 -= *a.offset((j + 1 as i32) as isize);
            j += 2 as i32
        }
        *a.offset((n - 1 as i32) as isize) = -xr
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_ooura_dfct(n: i32,
                                          a: *mut smpl_t,
                                          t: *mut smpl_t,
                                          ip: *mut i32,
                                          w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut l: i32 = 0;
    let mut m: i32 = 0;
    let mut mh: i32 = 0;
    let mut nw: i32 = 0;
    let mut nc: i32 = 0;
    let mut xr: smpl_t = 0.;
    let mut xi: smpl_t = 0.;
    let mut yr: smpl_t = 0.;
    let mut yi: smpl_t = 0.;
    nw = *ip.offset(0 as i32 as isize);
    if n > nw << 3 as i32 {
        nw = n >> 3 as i32;
        makewt(nw, ip, w);
    }
    nc = *ip.offset(1 as i32 as isize);
    if n > nc << 1 as i32 {
        nc = n >> 1 as i32;
        makect(nc, ip, w.offset(nw as isize));
    }
    m = n >> 1 as i32;
    yi = *a.offset(m as isize);
    xi = *a.offset(0 as i32 as isize) + *a.offset(n as isize);
    let ref mut fresh10 = *a.offset(0 as i32 as isize);
    *fresh10 -= *a.offset(n as isize);
    *t.offset(0 as i32 as isize) = xi - yi;
    *t.offset(m as isize) = xi + yi;
    if n > 2 as i32 {
        mh = m >> 1 as i32;
        j = 1 as i32;
        while j < mh {
            k = m - j;
            xr = *a.offset(j as isize) - *a.offset((n - j) as isize);
            xi = *a.offset(j as isize) + *a.offset((n - j) as isize);
            yr = *a.offset(k as isize) - *a.offset((n - k) as isize);
            yi = *a.offset(k as isize) + *a.offset((n - k) as isize);
            *a.offset(j as isize) = xr;
            *a.offset(k as isize) = yr;
            *t.offset(j as isize) = xi - yi;
            *t.offset(k as isize) = xi + yi;
            j += 1
        }
        *t.offset(mh as isize) =
            *a.offset(mh as isize) + *a.offset((n - mh) as isize);
        let ref mut fresh11 = *a.offset(mh as isize);
        *fresh11 -= *a.offset((n - mh) as isize);
        dctsub(m, a, nc, w.offset(nw as isize));
        if m > 4 as i32 {
            bitrv2(m, ip.offset(2 as i32 as isize), a);
            cftfsub(m, a, w);
            rftfsub(m, a, nc, w.offset(nw as isize));
        } else if m == 4 as i32 { cftfsub(m, a, w); }
        *a.offset((n - 1 as i32) as isize) =
            *a.offset(0 as i32 as isize) -
                *a.offset(1 as i32 as isize);
        *a.offset(1 as i32 as isize) =
            *a.offset(0 as i32 as isize) +
                *a.offset(1 as i32 as isize);
        j = m - 2 as i32;
        while j >= 2 as i32 {
            *a.offset((2 as i32 * j + 1 as i32) as isize) =
                *a.offset(j as isize) +
                    *a.offset((j + 1 as i32) as isize);
            *a.offset((2 as i32 * j - 1 as i32) as isize) =
                *a.offset(j as isize) -
                    *a.offset((j + 1 as i32) as isize);
            j -= 2 as i32
        }
        l = 2 as i32;
        m = mh;
        while m >= 2 as i32 {
            dctsub(m, t, nc, w.offset(nw as isize));
            if m > 4 as i32 {
                bitrv2(m, ip.offset(2 as i32 as isize), t);
                cftfsub(m, t, w);
                rftfsub(m, t, nc, w.offset(nw as isize));
            } else if m == 4 as i32 { cftfsub(m, t, w); }
            *a.offset((n - l) as isize) =
                *t.offset(0 as i32 as isize) -
                    *t.offset(1 as i32 as isize);
            *a.offset(l as isize) =
                *t.offset(0 as i32 as isize) +
                    *t.offset(1 as i32 as isize);
            k = 0 as i32;
            j = 2 as i32;
            while j < m {
                k += l << 2 as i32;
                *a.offset((k - l) as isize) =
                    *t.offset(j as isize) -
                        *t.offset((j + 1 as i32) as isize);
                *a.offset((k + l) as isize) =
                    *t.offset(j as isize) +
                        *t.offset((j + 1 as i32) as isize);
                j += 2 as i32
            }
            l <<= 1 as i32;
            mh = m >> 1 as i32;
            j = 0 as i32;
            while j < mh {
                k = m - j;
                *t.offset(j as isize) =
                    *t.offset((m + k) as isize) - *t.offset((m + j) as isize);
                *t.offset(k as isize) =
                    *t.offset((m + k) as isize) + *t.offset((m + j) as isize);
                j += 1
            }
            *t.offset(mh as isize) = *t.offset((m + mh) as isize);
            m = mh
        }
        *a.offset(l as isize) = *t.offset(0 as i32 as isize);
        *a.offset(n as isize) =
            *t.offset(2 as i32 as isize) -
                *t.offset(1 as i32 as isize);
        *a.offset(0 as i32 as isize) =
            *t.offset(2 as i32 as isize) +
                *t.offset(1 as i32 as isize)
    } else {
        *a.offset(1 as i32 as isize) =
            *a.offset(0 as i32 as isize);
        *a.offset(2 as i32 as isize) =
            *t.offset(0 as i32 as isize);
        *a.offset(0 as i32 as isize) =
            *t.offset(1 as i32 as isize)
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_ooura_dfst(n: i32,
                                          a: *mut smpl_t,
                                          t: *mut smpl_t,
                                          ip: *mut i32,
                                          w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut l: i32 = 0;
    let mut m: i32 = 0;
    let mut mh: i32 = 0;
    let mut nw: i32 = 0;
    let mut nc: i32 = 0;
    let mut xr: smpl_t = 0.;
    let mut xi: smpl_t = 0.;
    let mut yr: smpl_t = 0.;
    let mut yi: smpl_t = 0.;
    nw = *ip.offset(0 as i32 as isize);
    if n > nw << 3 as i32 {
        nw = n >> 3 as i32;
        makewt(nw, ip, w);
    }
    nc = *ip.offset(1 as i32 as isize);
    if n > nc << 1 as i32 {
        nc = n >> 1 as i32;
        makect(nc, ip, w.offset(nw as isize));
    }
    if n > 2 as i32 {
        m = n >> 1 as i32;
        mh = m >> 1 as i32;
        j = 1 as i32;
        while j < mh {
            k = m - j;
            xr = *a.offset(j as isize) + *a.offset((n - j) as isize);
            xi = *a.offset(j as isize) - *a.offset((n - j) as isize);
            yr = *a.offset(k as isize) + *a.offset((n - k) as isize);
            yi = *a.offset(k as isize) - *a.offset((n - k) as isize);
            *a.offset(j as isize) = xr;
            *a.offset(k as isize) = yr;
            *t.offset(j as isize) = xi + yi;
            *t.offset(k as isize) = xi - yi;
            j += 1
        }
        *t.offset(0 as i32 as isize) =
            *a.offset(mh as isize) - *a.offset((n - mh) as isize);
        let ref mut fresh12 = *a.offset(mh as isize);
        *fresh12 += *a.offset((n - mh) as isize);
        *a.offset(0 as i32 as isize) = *a.offset(m as isize);
        dstsub(m, a, nc, w.offset(nw as isize));
        if m > 4 as i32 {
            bitrv2(m, ip.offset(2 as i32 as isize), a);
            cftfsub(m, a, w);
            rftfsub(m, a, nc, w.offset(nw as isize));
        } else if m == 4 as i32 { cftfsub(m, a, w); }
        *a.offset((n - 1 as i32) as isize) =
            *a.offset(1 as i32 as isize) -
                *a.offset(0 as i32 as isize);
        *a.offset(1 as i32 as isize) =
            *a.offset(0 as i32 as isize) +
                *a.offset(1 as i32 as isize);
        j = m - 2 as i32;
        while j >= 2 as i32 {
            *a.offset((2 as i32 * j + 1 as i32) as isize) =
                *a.offset(j as isize) -
                    *a.offset((j + 1 as i32) as isize);
            *a.offset((2 as i32 * j - 1 as i32) as isize) =
                -*a.offset(j as isize) -
                    *a.offset((j + 1 as i32) as isize);
            j -= 2 as i32
        }
        l = 2 as i32;
        m = mh;
        while m >= 2 as i32 {
            dstsub(m, t, nc, w.offset(nw as isize));
            if m > 4 as i32 {
                bitrv2(m, ip.offset(2 as i32 as isize), t);
                cftfsub(m, t, w);
                rftfsub(m, t, nc, w.offset(nw as isize));
            } else if m == 4 as i32 { cftfsub(m, t, w); }
            *a.offset((n - l) as isize) =
                *t.offset(1 as i32 as isize) -
                    *t.offset(0 as i32 as isize);
            *a.offset(l as isize) =
                *t.offset(0 as i32 as isize) +
                    *t.offset(1 as i32 as isize);
            k = 0 as i32;
            j = 2 as i32;
            while j < m {
                k += l << 2 as i32;
                *a.offset((k - l) as isize) =
                    -*t.offset(j as isize) -
                        *t.offset((j + 1 as i32) as isize);
                *a.offset((k + l) as isize) =
                    *t.offset(j as isize) -
                        *t.offset((j + 1 as i32) as isize);
                j += 2 as i32
            }
            l <<= 1 as i32;
            mh = m >> 1 as i32;
            j = 1 as i32;
            while j < mh {
                k = m - j;
                *t.offset(j as isize) =
                    *t.offset((m + k) as isize) + *t.offset((m + j) as isize);
                *t.offset(k as isize) =
                    *t.offset((m + k) as isize) - *t.offset((m + j) as isize);
                j += 1
            }
            *t.offset(0 as i32 as isize) =
                *t.offset((m + mh) as isize);
            m = mh
        }
        *a.offset(l as isize) = *t.offset(0 as i32 as isize)
    }
    *a.offset(0 as i32 as isize) = 0 as i32 as smpl_t;
}
/* -------- initializing routines -------- */
unsafe extern "C" fn makewt(nw: i32, ip: *mut i32,
                            w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut nwh: i32 = 0;
    let mut delta: smpl_t = 0.;
    let mut x: smpl_t = 0.;
    let mut y: smpl_t = 0.;
    *ip.offset(0 as i32 as isize) = nw;
    *ip.offset(1 as i32 as isize) = 1 as i32;
    if nw > 2 as i32 {
        nwh = nw >> 1 as i32;
        delta = atanf(1.0f64 as f32) / nwh as f32;
        *w.offset(0 as i32 as isize) = 1 as i32 as smpl_t;
        *w.offset(1 as i32 as isize) = 0 as i32 as smpl_t;
        *w.offset(nwh as isize) = cosf(delta * nwh as f32);
        *w.offset((nwh + 1 as i32) as isize) =
            *w.offset(nwh as isize);
        if nwh > 2 as i32 {
            j = 2 as i32;
            while j < nwh {
                x = cosf(delta * j as f32);
                y = sinf(delta * j as f32);
                *w.offset(j as isize) = x;
                *w.offset((j + 1 as i32) as isize) = y;
                *w.offset((nw - j) as isize) = y;
                *w.offset((nw - j + 1 as i32) as isize) = x;
                j += 2 as i32
            }
            j = nwh - 2 as i32;
            while j >= 2 as i32 {
                x = *w.offset((2 as i32 * j) as isize);
                y =
                    *w.offset((2 as i32 * j + 1 as i32) as
                                  isize);
                *w.offset((nwh + j) as isize) = x;
                *w.offset((nwh + j + 1 as i32) as isize) = y;
                j -= 2 as i32
            }
            bitrv2(nw, ip.offset(2 as i32 as isize), w);
        }
    };
}
unsafe extern "C" fn makect(nc: i32, ip: *mut i32,
                            c: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut nch: i32 = 0;
    let mut delta: smpl_t = 0.;
    *ip.offset(1 as i32 as isize) = nc;
    if nc > 1 as i32 {
        nch = nc >> 1 as i32;
        delta = atanf(1.0f64 as f32) / nch as f32;
        *c.offset(0 as i32 as isize) =
            cosf(delta * nch as f32);
        *c.offset(nch as isize) =
            0.5f64 as smpl_t * *c.offset(0 as i32 as isize);
        j = 1 as i32;
        while j < nch {
            *c.offset(j as isize) =
                0.5f64 as smpl_t * cosf(delta * j as f32);
            *c.offset((nc - j) as isize) =
                0.5f64 as smpl_t * sinf(delta * j as f32);
            j += 1
        }
    };
}
/* -------- child routines -------- */
unsafe extern "C" fn bitrv2(n: i32, ip: *mut i32,
                            a: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut j1: i32 = 0;
    let mut k: i32 = 0;
    let mut k1: i32 = 0;
    let mut l: i32 = 0;
    let mut m: i32 = 0;
    let mut m2: i32 = 0;
    let mut xr: smpl_t = 0.;
    let mut xi: smpl_t = 0.;
    let mut yr: smpl_t = 0.;
    let mut yi: smpl_t = 0.;
    *ip.offset(0 as i32 as isize) = 0 as i32;
    l = n;
    m = 1 as i32;
    while (m << 3 as i32) < l {
        l >>= 1 as i32;
        j = 0 as i32;
        while j < m {
            *ip.offset((m + j) as isize) = *ip.offset(j as isize) + l;
            j += 1
        }
        m <<= 1 as i32
    }
    m2 = 2 as i32 * m;
    if m << 3 as i32 == l {
        k = 0 as i32;
        while k < m {
            j = 0 as i32;
            while j < k {
                j1 = 2 as i32 * j + *ip.offset(k as isize);
                k1 = 2 as i32 * k + *ip.offset(j as isize);
                xr = *a.offset(j1 as isize);
                xi = *a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = *a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 += 2 as i32 * m2;
                xr = *a.offset(j1 as isize);
                xi = *a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = *a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 -= m2;
                xr = *a.offset(j1 as isize);
                xi = *a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = *a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 += 2 as i32 * m2;
                xr = *a.offset(j1 as isize);
                xi = *a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = *a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j += 1
            }
            j1 = 2 as i32 * k + m2 + *ip.offset(k as isize);
            k1 = j1 + m2;
            xr = *a.offset(j1 as isize);
            xi = *a.offset((j1 + 1 as i32) as isize);
            yr = *a.offset(k1 as isize);
            yi = *a.offset((k1 + 1 as i32) as isize);
            *a.offset(j1 as isize) = yr;
            *a.offset((j1 + 1 as i32) as isize) = yi;
            *a.offset(k1 as isize) = xr;
            *a.offset((k1 + 1 as i32) as isize) = xi;
            k += 1
        }
    } else {
        k = 1 as i32;
        while k < m {
            j = 0 as i32;
            while j < k {
                j1 = 2 as i32 * j + *ip.offset(k as isize);
                k1 = 2 as i32 * k + *ip.offset(j as isize);
                xr = *a.offset(j1 as isize);
                xi = *a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = *a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 += m2;
                xr = *a.offset(j1 as isize);
                xi = *a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = *a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j += 1
            }
            k += 1
        }
    };
}
unsafe extern "C" fn bitrv2conj(n: i32, ip: *mut i32,
                                a: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut j1: i32 = 0;
    let mut k: i32 = 0;
    let mut k1: i32 = 0;
    let mut l: i32 = 0;
    let mut m: i32 = 0;
    let mut m2: i32 = 0;
    let mut xr: smpl_t = 0.;
    let mut xi: smpl_t = 0.;
    let mut yr: smpl_t = 0.;
    let mut yi: smpl_t = 0.;
    *ip.offset(0 as i32 as isize) = 0 as i32;
    l = n;
    m = 1 as i32;
    while (m << 3 as i32) < l {
        l >>= 1 as i32;
        j = 0 as i32;
        while j < m {
            *ip.offset((m + j) as isize) = *ip.offset(j as isize) + l;
            j += 1
        }
        m <<= 1 as i32
    }
    m2 = 2 as i32 * m;
    if m << 3 as i32 == l {
        k = 0 as i32;
        while k < m {
            j = 0 as i32;
            while j < k {
                j1 = 2 as i32 * j + *ip.offset(k as isize);
                k1 = 2 as i32 * k + *ip.offset(j as isize);
                xr = *a.offset(j1 as isize);
                xi = -*a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = -*a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 += 2 as i32 * m2;
                xr = *a.offset(j1 as isize);
                xi = -*a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = -*a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 -= m2;
                xr = *a.offset(j1 as isize);
                xi = -*a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = -*a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 += 2 as i32 * m2;
                xr = *a.offset(j1 as isize);
                xi = -*a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = -*a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j += 1
            }
            k1 = 2 as i32 * k + *ip.offset(k as isize);
            *a.offset((k1 + 1 as i32) as isize) =
                -*a.offset((k1 + 1 as i32) as isize);
            j1 = k1 + m2;
            k1 = j1 + m2;
            xr = *a.offset(j1 as isize);
            xi = -*a.offset((j1 + 1 as i32) as isize);
            yr = *a.offset(k1 as isize);
            yi = -*a.offset((k1 + 1 as i32) as isize);
            *a.offset(j1 as isize) = yr;
            *a.offset((j1 + 1 as i32) as isize) = yi;
            *a.offset(k1 as isize) = xr;
            *a.offset((k1 + 1 as i32) as isize) = xi;
            k1 += m2;
            *a.offset((k1 + 1 as i32) as isize) =
                -*a.offset((k1 + 1 as i32) as isize);
            k += 1
        }
    } else {
        *a.offset(1 as i32 as isize) =
            -*a.offset(1 as i32 as isize);
        *a.offset((m2 + 1 as i32) as isize) =
            -*a.offset((m2 + 1 as i32) as isize);
        k = 1 as i32;
        while k < m {
            j = 0 as i32;
            while j < k {
                j1 = 2 as i32 * j + *ip.offset(k as isize);
                k1 = 2 as i32 * k + *ip.offset(j as isize);
                xr = *a.offset(j1 as isize);
                xi = -*a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = -*a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j1 += m2;
                k1 += m2;
                xr = *a.offset(j1 as isize);
                xi = -*a.offset((j1 + 1 as i32) as isize);
                yr = *a.offset(k1 as isize);
                yi = -*a.offset((k1 + 1 as i32) as isize);
                *a.offset(j1 as isize) = yr;
                *a.offset((j1 + 1 as i32) as isize) = yi;
                *a.offset(k1 as isize) = xr;
                *a.offset((k1 + 1 as i32) as isize) = xi;
                j += 1
            }
            k1 = 2 as i32 * k + *ip.offset(k as isize);
            *a.offset((k1 + 1 as i32) as isize) =
                -*a.offset((k1 + 1 as i32) as isize);
            *a.offset((k1 + m2 + 1 as i32) as isize) =
                -*a.offset((k1 + m2 + 1 as i32) as isize);
            k += 1
        }
    };
}
unsafe extern "C" fn cftfsub(n: i32, a: *mut smpl_t,
                             w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut j1: i32 = 0;
    let mut j2: i32 = 0;
    let mut j3: i32 = 0;
    let mut l: i32 = 0;
    let mut x0r: smpl_t = 0.;
    let mut x0i: smpl_t = 0.;
    let mut x1r: smpl_t = 0.;
    let mut x1i: smpl_t = 0.;
    let mut x2r: smpl_t = 0.;
    let mut x2i: smpl_t = 0.;
    let mut x3r: smpl_t = 0.;
    let mut x3i: smpl_t = 0.;
    l = 2 as i32;
    if n >= 16 as i32 {
        cft1st(n, a, w);
        l = 16 as i32;
        while l << 3 as i32 <= n {
            cftmdl(n, l, a, w);
            l <<= 3 as i32
        }
    }
    if (l << 1 as i32) < n {
        j = 0 as i32;
        while j < l {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = *a.offset(j as isize) + *a.offset(j1 as isize);
            x0i =
                *a.offset((j + 1 as i32) as isize) +
                    *a.offset((j1 + 1 as i32) as isize);
            x1r = *a.offset(j as isize) - *a.offset(j1 as isize);
            x1i =
                *a.offset((j + 1 as i32) as isize) -
                    *a.offset((j1 + 1 as i32) as isize);
            x2r = *a.offset(j2 as isize) + *a.offset(j3 as isize);
            x2i =
                *a.offset((j2 + 1 as i32) as isize) +
                    *a.offset((j3 + 1 as i32) as isize);
            x3r = *a.offset(j2 as isize) - *a.offset(j3 as isize);
            x3i =
                *a.offset((j2 + 1 as i32) as isize) -
                    *a.offset((j3 + 1 as i32) as isize);
            *a.offset(j as isize) = x0r + x2r;
            *a.offset((j + 1 as i32) as isize) = x0i + x2i;
            *a.offset(j2 as isize) = x0r - x2r;
            *a.offset((j2 + 1 as i32) as isize) = x0i - x2i;
            *a.offset(j1 as isize) = x1r - x3i;
            *a.offset((j1 + 1 as i32) as isize) = x1i + x3r;
            *a.offset(j3 as isize) = x1r + x3i;
            *a.offset((j3 + 1 as i32) as isize) = x1i - x3r;
            j += 2 as i32
        }
    } else if l << 1 as i32 == n {
        j = 0 as i32;
        while j < l {
            j1 = j + l;
            x0r = *a.offset(j as isize) - *a.offset(j1 as isize);
            x0i =
                *a.offset((j + 1 as i32) as isize) -
                    *a.offset((j1 + 1 as i32) as isize);
            let ref mut fresh13 = *a.offset(j as isize);
            *fresh13 += *a.offset(j1 as isize);
            let ref mut fresh14 = *a.offset((j + 1 as i32) as isize);
            *fresh14 += *a.offset((j1 + 1 as i32) as isize);
            *a.offset(j1 as isize) = x0r;
            *a.offset((j1 + 1 as i32) as isize) = x0i;
            j += 2 as i32
        }
    };
}
unsafe extern "C" fn cftbsub(n: i32, a: *mut smpl_t,
                             w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut j1: i32 = 0;
    let mut j2: i32 = 0;
    let mut j3: i32 = 0;
    let mut j4: i32 = 0;
    let mut j5: i32 = 0;
    let mut j6: i32 = 0;
    let mut j7: i32 = 0;
    let mut l: i32 = 0;
    let mut wn4r: smpl_t = 0.;
    let mut x0r: smpl_t = 0.;
    let mut x0i: smpl_t = 0.;
    let mut x1r: smpl_t = 0.;
    let mut x1i: smpl_t = 0.;
    let mut x2r: smpl_t = 0.;
    let mut x2i: smpl_t = 0.;
    let mut x3r: smpl_t = 0.;
    let mut x3i: smpl_t = 0.;
    let mut y0r: smpl_t = 0.;
    let mut y0i: smpl_t = 0.;
    let mut y1r: smpl_t = 0.;
    let mut y1i: smpl_t = 0.;
    let mut y2r: smpl_t = 0.;
    let mut y2i: smpl_t = 0.;
    let mut y3r: smpl_t = 0.;
    let mut y3i: smpl_t = 0.;
    let mut y4r: smpl_t = 0.;
    let mut y4i: smpl_t = 0.;
    let mut y5r: smpl_t = 0.;
    let mut y5i: smpl_t = 0.;
    let mut y6r: smpl_t = 0.;
    let mut y6i: smpl_t = 0.;
    let mut y7r: smpl_t = 0.;
    let mut y7i: smpl_t = 0.;
    l = 2 as i32;
    if n > 16 as i32 {
        cft1st(n, a, w);
        l = 16 as i32;
        while (l << 3 as i32) < n {
            cftmdl(n, l, a, w);
            l <<= 3 as i32
        }
    }
    if (l << 2 as i32) < n {
        wn4r = *w.offset(2 as i32 as isize);
        j = 0 as i32;
        while j < l {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            j4 = j3 + l;
            j5 = j4 + l;
            j6 = j5 + l;
            j7 = j6 + l;
            x0r = *a.offset(j as isize) + *a.offset(j1 as isize);
            x0i =
                -*a.offset((j + 1 as i32) as isize) -
                    *a.offset((j1 + 1 as i32) as isize);
            x1r = *a.offset(j as isize) - *a.offset(j1 as isize);
            x1i =
                -*a.offset((j + 1 as i32) as isize) +
                    *a.offset((j1 + 1 as i32) as isize);
            x2r = *a.offset(j2 as isize) + *a.offset(j3 as isize);
            x2i =
                *a.offset((j2 + 1 as i32) as isize) +
                    *a.offset((j3 + 1 as i32) as isize);
            x3r = *a.offset(j2 as isize) - *a.offset(j3 as isize);
            x3i =
                *a.offset((j2 + 1 as i32) as isize) -
                    *a.offset((j3 + 1 as i32) as isize);
            y0r = x0r + x2r;
            y0i = x0i - x2i;
            y2r = x0r - x2r;
            y2i = x0i + x2i;
            y1r = x1r - x3i;
            y1i = x1i - x3r;
            y3r = x1r + x3i;
            y3i = x1i + x3r;
            x0r = *a.offset(j4 as isize) + *a.offset(j5 as isize);
            x0i =
                *a.offset((j4 + 1 as i32) as isize) +
                    *a.offset((j5 + 1 as i32) as isize);
            x1r = *a.offset(j4 as isize) - *a.offset(j5 as isize);
            x1i =
                *a.offset((j4 + 1 as i32) as isize) -
                    *a.offset((j5 + 1 as i32) as isize);
            x2r = *a.offset(j6 as isize) + *a.offset(j7 as isize);
            x2i =
                *a.offset((j6 + 1 as i32) as isize) +
                    *a.offset((j7 + 1 as i32) as isize);
            x3r = *a.offset(j6 as isize) - *a.offset(j7 as isize);
            x3i =
                *a.offset((j6 + 1 as i32) as isize) -
                    *a.offset((j7 + 1 as i32) as isize);
            y4r = x0r + x2r;
            y4i = x0i + x2i;
            y6r = x0r - x2r;
            y6i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            x2r = x1r + x3i;
            x2i = x1i - x3r;
            y5r = wn4r * (x0r - x0i);
            y5i = wn4r * (x0r + x0i);
            y7r = wn4r * (x2r - x2i);
            y7i = wn4r * (x2r + x2i);
            *a.offset(j1 as isize) = y1r + y5r;
            *a.offset((j1 + 1 as i32) as isize) = y1i - y5i;
            *a.offset(j5 as isize) = y1r - y5r;
            *a.offset((j5 + 1 as i32) as isize) = y1i + y5i;
            *a.offset(j3 as isize) = y3r - y7i;
            *a.offset((j3 + 1 as i32) as isize) = y3i - y7r;
            *a.offset(j7 as isize) = y3r + y7i;
            *a.offset((j7 + 1 as i32) as isize) = y3i + y7r;
            *a.offset(j as isize) = y0r + y4r;
            *a.offset((j + 1 as i32) as isize) = y0i - y4i;
            *a.offset(j4 as isize) = y0r - y4r;
            *a.offset((j4 + 1 as i32) as isize) = y0i + y4i;
            *a.offset(j2 as isize) = y2r - y6i;
            *a.offset((j2 + 1 as i32) as isize) = y2i - y6r;
            *a.offset(j6 as isize) = y2r + y6i;
            *a.offset((j6 + 1 as i32) as isize) = y2i + y6r;
            j += 2 as i32
        }
    } else if l << 2 as i32 == n {
        j = 0 as i32;
        while j < l {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = *a.offset(j as isize) + *a.offset(j1 as isize);
            x0i =
                -*a.offset((j + 1 as i32) as isize) -
                    *a.offset((j1 + 1 as i32) as isize);
            x1r = *a.offset(j as isize) - *a.offset(j1 as isize);
            x1i =
                -*a.offset((j + 1 as i32) as isize) +
                    *a.offset((j1 + 1 as i32) as isize);
            x2r = *a.offset(j2 as isize) + *a.offset(j3 as isize);
            x2i =
                *a.offset((j2 + 1 as i32) as isize) +
                    *a.offset((j3 + 1 as i32) as isize);
            x3r = *a.offset(j2 as isize) - *a.offset(j3 as isize);
            x3i =
                *a.offset((j2 + 1 as i32) as isize) -
                    *a.offset((j3 + 1 as i32) as isize);
            *a.offset(j as isize) = x0r + x2r;
            *a.offset((j + 1 as i32) as isize) = x0i - x2i;
            *a.offset(j2 as isize) = x0r - x2r;
            *a.offset((j2 + 1 as i32) as isize) = x0i + x2i;
            *a.offset(j1 as isize) = x1r - x3i;
            *a.offset((j1 + 1 as i32) as isize) = x1i - x3r;
            *a.offset(j3 as isize) = x1r + x3i;
            *a.offset((j3 + 1 as i32) as isize) = x1i + x3r;
            j += 2 as i32
        }
    } else {
        j = 0 as i32;
        while j < l {
            j1 = j + l;
            x0r = *a.offset(j as isize) - *a.offset(j1 as isize);
            x0i =
                -*a.offset((j + 1 as i32) as isize) +
                    *a.offset((j1 + 1 as i32) as isize);
            let ref mut fresh15 = *a.offset(j as isize);
            *fresh15 += *a.offset(j1 as isize);
            *a.offset((j + 1 as i32) as isize) =
                -*a.offset((j + 1 as i32) as isize) -
                    *a.offset((j1 + 1 as i32) as isize);
            *a.offset(j1 as isize) = x0r;
            *a.offset((j1 + 1 as i32) as isize) = x0i;
            j += 2 as i32
        }
    };
}
unsafe extern "C" fn cft1st(n: i32, a: *mut smpl_t,
                            w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k1: i32 = 0;
    let mut wn4r: smpl_t = 0.;
    let mut wtmp: smpl_t = 0.;
    let mut wk1r: smpl_t = 0.;
    let mut wk1i: smpl_t = 0.;
    let mut wk2r: smpl_t = 0.;
    let mut wk2i: smpl_t = 0.;
    let mut wk3r: smpl_t = 0.;
    let mut wk3i: smpl_t = 0.;
    let mut wk4r: smpl_t = 0.;
    let mut wk4i: smpl_t = 0.;
    let mut wk5r: smpl_t = 0.;
    let mut wk5i: smpl_t = 0.;
    let mut wk6r: smpl_t = 0.;
    let mut wk6i: smpl_t = 0.;
    let mut wk7r: smpl_t = 0.;
    let mut wk7i: smpl_t = 0.;
    let mut x0r: smpl_t = 0.;
    let mut x0i: smpl_t = 0.;
    let mut x1r: smpl_t = 0.;
    let mut x1i: smpl_t = 0.;
    let mut x2r: smpl_t = 0.;
    let mut x2i: smpl_t = 0.;
    let mut x3r: smpl_t = 0.;
    let mut x3i: smpl_t = 0.;
    let mut y0r: smpl_t = 0.;
    let mut y0i: smpl_t = 0.;
    let mut y1r: smpl_t = 0.;
    let mut y1i: smpl_t = 0.;
    let mut y2r: smpl_t = 0.;
    let mut y2i: smpl_t = 0.;
    let mut y3r: smpl_t = 0.;
    let mut y3i: smpl_t = 0.;
    let mut y4r: smpl_t = 0.;
    let mut y4i: smpl_t = 0.;
    let mut y5r: smpl_t = 0.;
    let mut y5i: smpl_t = 0.;
    let mut y6r: smpl_t = 0.;
    let mut y6i: smpl_t = 0.;
    let mut y7r: smpl_t = 0.;
    let mut y7i: smpl_t = 0.;
    wn4r = *w.offset(2 as i32 as isize);
    x0r =
        *a.offset(0 as i32 as isize) +
            *a.offset(2 as i32 as isize);
    x0i =
        *a.offset(1 as i32 as isize) +
            *a.offset(3 as i32 as isize);
    x1r =
        *a.offset(0 as i32 as isize) -
            *a.offset(2 as i32 as isize);
    x1i =
        *a.offset(1 as i32 as isize) -
            *a.offset(3 as i32 as isize);
    x2r =
        *a.offset(4 as i32 as isize) +
            *a.offset(6 as i32 as isize);
    x2i =
        *a.offset(5 as i32 as isize) +
            *a.offset(7 as i32 as isize);
    x3r =
        *a.offset(4 as i32 as isize) -
            *a.offset(6 as i32 as isize);
    x3i =
        *a.offset(5 as i32 as isize) -
            *a.offset(7 as i32 as isize);
    y0r = x0r + x2r;
    y0i = x0i + x2i;
    y2r = x0r - x2r;
    y2i = x0i - x2i;
    y1r = x1r - x3i;
    y1i = x1i + x3r;
    y3r = x1r + x3i;
    y3i = x1i - x3r;
    x0r =
        *a.offset(8 as i32 as isize) +
            *a.offset(10 as i32 as isize);
    x0i =
        *a.offset(9 as i32 as isize) +
            *a.offset(11 as i32 as isize);
    x1r =
        *a.offset(8 as i32 as isize) -
            *a.offset(10 as i32 as isize);
    x1i =
        *a.offset(9 as i32 as isize) -
            *a.offset(11 as i32 as isize);
    x2r =
        *a.offset(12 as i32 as isize) +
            *a.offset(14 as i32 as isize);
    x2i =
        *a.offset(13 as i32 as isize) +
            *a.offset(15 as i32 as isize);
    x3r =
        *a.offset(12 as i32 as isize) -
            *a.offset(14 as i32 as isize);
    x3i =
        *a.offset(13 as i32 as isize) -
            *a.offset(15 as i32 as isize);
    y4r = x0r + x2r;
    y4i = x0i + x2i;
    y6r = x0r - x2r;
    y6i = x0i - x2i;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    x2r = x1r + x3i;
    x2i = x1i - x3r;
    y5r = wn4r * (x0r - x0i);
    y5i = wn4r * (x0r + x0i);
    y7r = wn4r * (x2r - x2i);
    y7i = wn4r * (x2r + x2i);
    *a.offset(2 as i32 as isize) = y1r + y5r;
    *a.offset(3 as i32 as isize) = y1i + y5i;
    *a.offset(10 as i32 as isize) = y1r - y5r;
    *a.offset(11 as i32 as isize) = y1i - y5i;
    *a.offset(6 as i32 as isize) = y3r - y7i;
    *a.offset(7 as i32 as isize) = y3i + y7r;
    *a.offset(14 as i32 as isize) = y3r + y7i;
    *a.offset(15 as i32 as isize) = y3i - y7r;
    *a.offset(0 as i32 as isize) = y0r + y4r;
    *a.offset(1 as i32 as isize) = y0i + y4i;
    *a.offset(8 as i32 as isize) = y0r - y4r;
    *a.offset(9 as i32 as isize) = y0i - y4i;
    *a.offset(4 as i32 as isize) = y2r - y6i;
    *a.offset(5 as i32 as isize) = y2i + y6r;
    *a.offset(12 as i32 as isize) = y2r + y6i;
    *a.offset(13 as i32 as isize) = y2i - y6r;
    if n > 16 as i32 {
        wk1r = *w.offset(4 as i32 as isize);
        wk1i = *w.offset(5 as i32 as isize);
        x0r =
            *a.offset(16 as i32 as isize) +
                *a.offset(18 as i32 as isize);
        x0i =
            *a.offset(17 as i32 as isize) +
                *a.offset(19 as i32 as isize);
        x1r =
            *a.offset(16 as i32 as isize) -
                *a.offset(18 as i32 as isize);
        x1i =
            *a.offset(17 as i32 as isize) -
                *a.offset(19 as i32 as isize);
        x2r =
            *a.offset(20 as i32 as isize) +
                *a.offset(22 as i32 as isize);
        x2i =
            *a.offset(21 as i32 as isize) +
                *a.offset(23 as i32 as isize);
        x3r =
            *a.offset(20 as i32 as isize) -
                *a.offset(22 as i32 as isize);
        x3i =
            *a.offset(21 as i32 as isize) -
                *a.offset(23 as i32 as isize);
        y0r = x0r + x2r;
        y0i = x0i + x2i;
        y2r = x0r - x2r;
        y2i = x0i - x2i;
        y1r = x1r - x3i;
        y1i = x1i + x3r;
        y3r = x1r + x3i;
        y3i = x1i - x3r;
        x0r =
            *a.offset(24 as i32 as isize) +
                *a.offset(26 as i32 as isize);
        x0i =
            *a.offset(25 as i32 as isize) +
                *a.offset(27 as i32 as isize);
        x1r =
            *a.offset(24 as i32 as isize) -
                *a.offset(26 as i32 as isize);
        x1i =
            *a.offset(25 as i32 as isize) -
                *a.offset(27 as i32 as isize);
        x2r =
            *a.offset(28 as i32 as isize) +
                *a.offset(30 as i32 as isize);
        x2i =
            *a.offset(29 as i32 as isize) +
                *a.offset(31 as i32 as isize);
        x3r =
            *a.offset(28 as i32 as isize) -
                *a.offset(30 as i32 as isize);
        x3i =
            *a.offset(29 as i32 as isize) -
                *a.offset(31 as i32 as isize);
        y4r = x0r + x2r;
        y4i = x0i + x2i;
        y6r = x0r - x2r;
        y6i = x0i - x2i;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        x2r = x1r + x3i;
        x2i = x3r - x1i;
        y5r = wk1i * x0r - wk1r * x0i;
        y5i = wk1i * x0i + wk1r * x0r;
        y7r = wk1r * x2r + wk1i * x2i;
        y7i = wk1r * x2i - wk1i * x2r;
        x0r = wk1r * y1r - wk1i * y1i;
        x0i = wk1r * y1i + wk1i * y1r;
        *a.offset(18 as i32 as isize) = x0r + y5r;
        *a.offset(19 as i32 as isize) = x0i + y5i;
        *a.offset(26 as i32 as isize) = y5i - x0i;
        *a.offset(27 as i32 as isize) = x0r - y5r;
        x0r = wk1i * y3r - wk1r * y3i;
        x0i = wk1i * y3i + wk1r * y3r;
        *a.offset(22 as i32 as isize) = x0r - y7r;
        *a.offset(23 as i32 as isize) = x0i + y7i;
        *a.offset(30 as i32 as isize) = y7i - x0i;
        *a.offset(31 as i32 as isize) = x0r + y7r;
        *a.offset(16 as i32 as isize) = y0r + y4r;
        *a.offset(17 as i32 as isize) = y0i + y4i;
        *a.offset(24 as i32 as isize) = y4i - y0i;
        *a.offset(25 as i32 as isize) = y0r - y4r;
        x0r = y2r - y6i;
        x0i = y2i + y6r;
        *a.offset(20 as i32 as isize) = wn4r * (x0r - x0i);
        *a.offset(21 as i32 as isize) = wn4r * (x0i + x0r);
        x0r = y6r - y2i;
        x0i = y2r + y6i;
        *a.offset(28 as i32 as isize) = wn4r * (x0r - x0i);
        *a.offset(29 as i32 as isize) = wn4r * (x0i + x0r);
        k1 = 4 as i32;
        j = 32 as i32;
        while j < n {
            k1 += 4 as i32;
            wk1r = *w.offset(k1 as isize);
            wk1i = *w.offset((k1 + 1 as i32) as isize);
            wk2r = *w.offset((k1 + 2 as i32) as isize);
            wk2i = *w.offset((k1 + 3 as i32) as isize);
            wtmp = 2 as i32 as f32 * wk2i;
            wk3r = wk1r - wtmp * wk1i;
            wk3i = wtmp * wk1r - wk1i;
            wk4r = 1 as i32 as f32 - wtmp * wk2i;
            wk4i = wtmp * wk2r;
            wtmp = 2 as i32 as f32 * wk4i;
            wk5r = wk3r - wtmp * wk1i;
            wk5i = wtmp * wk1r - wk3i;
            wk6r = wk2r - wtmp * wk2i;
            wk6i = wtmp * wk2r - wk2i;
            wk7r = wk1r - wtmp * wk3i;
            wk7i = wtmp * wk3r - wk1i;
            x0r =
                *a.offset(j as isize) +
                    *a.offset((j + 2 as i32) as isize);
            x0i =
                *a.offset((j + 1 as i32) as isize) +
                    *a.offset((j + 3 as i32) as isize);
            x1r =
                *a.offset(j as isize) -
                    *a.offset((j + 2 as i32) as isize);
            x1i =
                *a.offset((j + 1 as i32) as isize) -
                    *a.offset((j + 3 as i32) as isize);
            x2r =
                *a.offset((j + 4 as i32) as isize) +
                    *a.offset((j + 6 as i32) as isize);
            x2i =
                *a.offset((j + 5 as i32) as isize) +
                    *a.offset((j + 7 as i32) as isize);
            x3r =
                *a.offset((j + 4 as i32) as isize) -
                    *a.offset((j + 6 as i32) as isize);
            x3i =
                *a.offset((j + 5 as i32) as isize) -
                    *a.offset((j + 7 as i32) as isize);
            y0r = x0r + x2r;
            y0i = x0i + x2i;
            y2r = x0r - x2r;
            y2i = x0i - x2i;
            y1r = x1r - x3i;
            y1i = x1i + x3r;
            y3r = x1r + x3i;
            y3i = x1i - x3r;
            x0r =
                *a.offset((j + 8 as i32) as isize) +
                    *a.offset((j + 10 as i32) as isize);
            x0i =
                *a.offset((j + 9 as i32) as isize) +
                    *a.offset((j + 11 as i32) as isize);
            x1r =
                *a.offset((j + 8 as i32) as isize) -
                    *a.offset((j + 10 as i32) as isize);
            x1i =
                *a.offset((j + 9 as i32) as isize) -
                    *a.offset((j + 11 as i32) as isize);
            x2r =
                *a.offset((j + 12 as i32) as isize) +
                    *a.offset((j + 14 as i32) as isize);
            x2i =
                *a.offset((j + 13 as i32) as isize) +
                    *a.offset((j + 15 as i32) as isize);
            x3r =
                *a.offset((j + 12 as i32) as isize) -
                    *a.offset((j + 14 as i32) as isize);
            x3i =
                *a.offset((j + 13 as i32) as isize) -
                    *a.offset((j + 15 as i32) as isize);
            y4r = x0r + x2r;
            y4i = x0i + x2i;
            y6r = x0r - x2r;
            y6i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            x2r = x1r + x3i;
            x2i = x1i - x3r;
            y5r = wn4r * (x0r - x0i);
            y5i = wn4r * (x0r + x0i);
            y7r = wn4r * (x2r - x2i);
            y7i = wn4r * (x2r + x2i);
            x0r = y1r + y5r;
            x0i = y1i + y5i;
            *a.offset((j + 2 as i32) as isize) =
                wk1r * x0r - wk1i * x0i;
            *a.offset((j + 3 as i32) as isize) =
                wk1r * x0i + wk1i * x0r;
            x0r = y1r - y5r;
            x0i = y1i - y5i;
            *a.offset((j + 10 as i32) as isize) =
                wk5r * x0r - wk5i * x0i;
            *a.offset((j + 11 as i32) as isize) =
                wk5r * x0i + wk5i * x0r;
            x0r = y3r - y7i;
            x0i = y3i + y7r;
            *a.offset((j + 6 as i32) as isize) =
                wk3r * x0r - wk3i * x0i;
            *a.offset((j + 7 as i32) as isize) =
                wk3r * x0i + wk3i * x0r;
            x0r = y3r + y7i;
            x0i = y3i - y7r;
            *a.offset((j + 14 as i32) as isize) =
                wk7r * x0r - wk7i * x0i;
            *a.offset((j + 15 as i32) as isize) =
                wk7r * x0i + wk7i * x0r;
            *a.offset(j as isize) = y0r + y4r;
            *a.offset((j + 1 as i32) as isize) = y0i + y4i;
            x0r = y0r - y4r;
            x0i = y0i - y4i;
            *a.offset((j + 8 as i32) as isize) =
                wk4r * x0r - wk4i * x0i;
            *a.offset((j + 9 as i32) as isize) =
                wk4r * x0i + wk4i * x0r;
            x0r = y2r - y6i;
            x0i = y2i + y6r;
            *a.offset((j + 4 as i32) as isize) =
                wk2r * x0r - wk2i * x0i;
            *a.offset((j + 5 as i32) as isize) =
                wk2r * x0i + wk2i * x0r;
            x0r = y2r + y6i;
            x0i = y2i - y6r;
            *a.offset((j + 12 as i32) as isize) =
                wk6r * x0r - wk6i * x0i;
            *a.offset((j + 13 as i32) as isize) =
                wk6r * x0i + wk6i * x0r;
            j += 16 as i32
        }
    };
}
unsafe extern "C" fn cftmdl(n: i32, l: i32,
                            a: *mut smpl_t, w: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut j1: i32 = 0;
    let mut j2: i32 = 0;
    let mut j3: i32 = 0;
    let mut j4: i32 = 0;
    let mut j5: i32 = 0;
    let mut j6: i32 = 0;
    let mut j7: i32 = 0;
    let mut k: i32 = 0;
    let mut k1: i32 = 0;
    let mut m: i32 = 0;
    let mut wn4r: smpl_t = 0.;
    let mut wtmp: smpl_t = 0.;
    let mut wk1r: smpl_t = 0.;
    let mut wk1i: smpl_t = 0.;
    let mut wk2r: smpl_t = 0.;
    let mut wk2i: smpl_t = 0.;
    let mut wk3r: smpl_t = 0.;
    let mut wk3i: smpl_t = 0.;
    let mut wk4r: smpl_t = 0.;
    let mut wk4i: smpl_t = 0.;
    let mut wk5r: smpl_t = 0.;
    let mut wk5i: smpl_t = 0.;
    let mut wk6r: smpl_t = 0.;
    let mut wk6i: smpl_t = 0.;
    let mut wk7r: smpl_t = 0.;
    let mut wk7i: smpl_t = 0.;
    let mut x0r: smpl_t = 0.;
    let mut x0i: smpl_t = 0.;
    let mut x1r: smpl_t = 0.;
    let mut x1i: smpl_t = 0.;
    let mut x2r: smpl_t = 0.;
    let mut x2i: smpl_t = 0.;
    let mut x3r: smpl_t = 0.;
    let mut x3i: smpl_t = 0.;
    let mut y0r: smpl_t = 0.;
    let mut y0i: smpl_t = 0.;
    let mut y1r: smpl_t = 0.;
    let mut y1i: smpl_t = 0.;
    let mut y2r: smpl_t = 0.;
    let mut y2i: smpl_t = 0.;
    let mut y3r: smpl_t = 0.;
    let mut y3i: smpl_t = 0.;
    let mut y4r: smpl_t = 0.;
    let mut y4i: smpl_t = 0.;
    let mut y5r: smpl_t = 0.;
    let mut y5i: smpl_t = 0.;
    let mut y6r: smpl_t = 0.;
    let mut y6i: smpl_t = 0.;
    let mut y7r: smpl_t = 0.;
    let mut y7i: smpl_t = 0.;
    m = l << 3 as i32;
    wn4r = *w.offset(2 as i32 as isize);
    j = 0 as i32;
    while j < l {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        j4 = j3 + l;
        j5 = j4 + l;
        j6 = j5 + l;
        j7 = j6 + l;
        x0r = *a.offset(j as isize) + *a.offset(j1 as isize);
        x0i =
            *a.offset((j + 1 as i32) as isize) +
                *a.offset((j1 + 1 as i32) as isize);
        x1r = *a.offset(j as isize) - *a.offset(j1 as isize);
        x1i =
            *a.offset((j + 1 as i32) as isize) -
                *a.offset((j1 + 1 as i32) as isize);
        x2r = *a.offset(j2 as isize) + *a.offset(j3 as isize);
        x2i =
            *a.offset((j2 + 1 as i32) as isize) +
                *a.offset((j3 + 1 as i32) as isize);
        x3r = *a.offset(j2 as isize) - *a.offset(j3 as isize);
        x3i =
            *a.offset((j2 + 1 as i32) as isize) -
                *a.offset((j3 + 1 as i32) as isize);
        y0r = x0r + x2r;
        y0i = x0i + x2i;
        y2r = x0r - x2r;
        y2i = x0i - x2i;
        y1r = x1r - x3i;
        y1i = x1i + x3r;
        y3r = x1r + x3i;
        y3i = x1i - x3r;
        x0r = *a.offset(j4 as isize) + *a.offset(j5 as isize);
        x0i =
            *a.offset((j4 + 1 as i32) as isize) +
                *a.offset((j5 + 1 as i32) as isize);
        x1r = *a.offset(j4 as isize) - *a.offset(j5 as isize);
        x1i =
            *a.offset((j4 + 1 as i32) as isize) -
                *a.offset((j5 + 1 as i32) as isize);
        x2r = *a.offset(j6 as isize) + *a.offset(j7 as isize);
        x2i =
            *a.offset((j6 + 1 as i32) as isize) +
                *a.offset((j7 + 1 as i32) as isize);
        x3r = *a.offset(j6 as isize) - *a.offset(j7 as isize);
        x3i =
            *a.offset((j6 + 1 as i32) as isize) -
                *a.offset((j7 + 1 as i32) as isize);
        y4r = x0r + x2r;
        y4i = x0i + x2i;
        y6r = x0r - x2r;
        y6i = x0i - x2i;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        x2r = x1r + x3i;
        x2i = x1i - x3r;
        y5r = wn4r * (x0r - x0i);
        y5i = wn4r * (x0r + x0i);
        y7r = wn4r * (x2r - x2i);
        y7i = wn4r * (x2r + x2i);
        *a.offset(j1 as isize) = y1r + y5r;
        *a.offset((j1 + 1 as i32) as isize) = y1i + y5i;
        *a.offset(j5 as isize) = y1r - y5r;
        *a.offset((j5 + 1 as i32) as isize) = y1i - y5i;
        *a.offset(j3 as isize) = y3r - y7i;
        *a.offset((j3 + 1 as i32) as isize) = y3i + y7r;
        *a.offset(j7 as isize) = y3r + y7i;
        *a.offset((j7 + 1 as i32) as isize) = y3i - y7r;
        *a.offset(j as isize) = y0r + y4r;
        *a.offset((j + 1 as i32) as isize) = y0i + y4i;
        *a.offset(j4 as isize) = y0r - y4r;
        *a.offset((j4 + 1 as i32) as isize) = y0i - y4i;
        *a.offset(j2 as isize) = y2r - y6i;
        *a.offset((j2 + 1 as i32) as isize) = y2i + y6r;
        *a.offset(j6 as isize) = y2r + y6i;
        *a.offset((j6 + 1 as i32) as isize) = y2i - y6r;
        j += 2 as i32
    }
    if m < n {
        wk1r = *w.offset(4 as i32 as isize);
        wk1i = *w.offset(5 as i32 as isize);
        j = m;
        while j < l + m {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            j4 = j3 + l;
            j5 = j4 + l;
            j6 = j5 + l;
            j7 = j6 + l;
            x0r = *a.offset(j as isize) + *a.offset(j1 as isize);
            x0i =
                *a.offset((j + 1 as i32) as isize) +
                    *a.offset((j1 + 1 as i32) as isize);
            x1r = *a.offset(j as isize) - *a.offset(j1 as isize);
            x1i =
                *a.offset((j + 1 as i32) as isize) -
                    *a.offset((j1 + 1 as i32) as isize);
            x2r = *a.offset(j2 as isize) + *a.offset(j3 as isize);
            x2i =
                *a.offset((j2 + 1 as i32) as isize) +
                    *a.offset((j3 + 1 as i32) as isize);
            x3r = *a.offset(j2 as isize) - *a.offset(j3 as isize);
            x3i =
                *a.offset((j2 + 1 as i32) as isize) -
                    *a.offset((j3 + 1 as i32) as isize);
            y0r = x0r + x2r;
            y0i = x0i + x2i;
            y2r = x0r - x2r;
            y2i = x0i - x2i;
            y1r = x1r - x3i;
            y1i = x1i + x3r;
            y3r = x1r + x3i;
            y3i = x1i - x3r;
            x0r = *a.offset(j4 as isize) + *a.offset(j5 as isize);
            x0i =
                *a.offset((j4 + 1 as i32) as isize) +
                    *a.offset((j5 + 1 as i32) as isize);
            x1r = *a.offset(j4 as isize) - *a.offset(j5 as isize);
            x1i =
                *a.offset((j4 + 1 as i32) as isize) -
                    *a.offset((j5 + 1 as i32) as isize);
            x2r = *a.offset(j6 as isize) + *a.offset(j7 as isize);
            x2i =
                *a.offset((j6 + 1 as i32) as isize) +
                    *a.offset((j7 + 1 as i32) as isize);
            x3r = *a.offset(j6 as isize) - *a.offset(j7 as isize);
            x3i =
                *a.offset((j6 + 1 as i32) as isize) -
                    *a.offset((j7 + 1 as i32) as isize);
            y4r = x0r + x2r;
            y4i = x0i + x2i;
            y6r = x0r - x2r;
            y6i = x0i - x2i;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            x2r = x1r + x3i;
            x2i = x3r - x1i;
            y5r = wk1i * x0r - wk1r * x0i;
            y5i = wk1i * x0i + wk1r * x0r;
            y7r = wk1r * x2r + wk1i * x2i;
            y7i = wk1r * x2i - wk1i * x2r;
            x0r = wk1r * y1r - wk1i * y1i;
            x0i = wk1r * y1i + wk1i * y1r;
            *a.offset(j1 as isize) = x0r + y5r;
            *a.offset((j1 + 1 as i32) as isize) = x0i + y5i;
            *a.offset(j5 as isize) = y5i - x0i;
            *a.offset((j5 + 1 as i32) as isize) = x0r - y5r;
            x0r = wk1i * y3r - wk1r * y3i;
            x0i = wk1i * y3i + wk1r * y3r;
            *a.offset(j3 as isize) = x0r - y7r;
            *a.offset((j3 + 1 as i32) as isize) = x0i + y7i;
            *a.offset(j7 as isize) = y7i - x0i;
            *a.offset((j7 + 1 as i32) as isize) = x0r + y7r;
            *a.offset(j as isize) = y0r + y4r;
            *a.offset((j + 1 as i32) as isize) = y0i + y4i;
            *a.offset(j4 as isize) = y4i - y0i;
            *a.offset((j4 + 1 as i32) as isize) = y0r - y4r;
            x0r = y2r - y6i;
            x0i = y2i + y6r;
            *a.offset(j2 as isize) = wn4r * (x0r - x0i);
            *a.offset((j2 + 1 as i32) as isize) = wn4r * (x0i + x0r);
            x0r = y6r - y2i;
            x0i = y2r + y6i;
            *a.offset(j6 as isize) = wn4r * (x0r - x0i);
            *a.offset((j6 + 1 as i32) as isize) = wn4r * (x0i + x0r);
            j += 2 as i32
        }
        k1 = 4 as i32;
        k = 2 as i32 * m;
        while k < n {
            k1 += 4 as i32;
            wk1r = *w.offset(k1 as isize);
            wk1i = *w.offset((k1 + 1 as i32) as isize);
            wk2r = *w.offset((k1 + 2 as i32) as isize);
            wk2i = *w.offset((k1 + 3 as i32) as isize);
            wtmp = 2 as i32 as f32 * wk2i;
            wk3r = wk1r - wtmp * wk1i;
            wk3i = wtmp * wk1r - wk1i;
            wk4r = 1 as i32 as f32 - wtmp * wk2i;
            wk4i = wtmp * wk2r;
            wtmp = 2 as i32 as f32 * wk4i;
            wk5r = wk3r - wtmp * wk1i;
            wk5i = wtmp * wk1r - wk3i;
            wk6r = wk2r - wtmp * wk2i;
            wk6i = wtmp * wk2r - wk2i;
            wk7r = wk1r - wtmp * wk3i;
            wk7i = wtmp * wk3r - wk1i;
            j = k;
            while j < l + k {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                j4 = j3 + l;
                j5 = j4 + l;
                j6 = j5 + l;
                j7 = j6 + l;
                x0r = *a.offset(j as isize) + *a.offset(j1 as isize);
                x0i =
                    *a.offset((j + 1 as i32) as isize) +
                        *a.offset((j1 + 1 as i32) as isize);
                x1r = *a.offset(j as isize) - *a.offset(j1 as isize);
                x1i =
                    *a.offset((j + 1 as i32) as isize) -
                        *a.offset((j1 + 1 as i32) as isize);
                x2r = *a.offset(j2 as isize) + *a.offset(j3 as isize);
                x2i =
                    *a.offset((j2 + 1 as i32) as isize) +
                        *a.offset((j3 + 1 as i32) as isize);
                x3r = *a.offset(j2 as isize) - *a.offset(j3 as isize);
                x3i =
                    *a.offset((j2 + 1 as i32) as isize) -
                        *a.offset((j3 + 1 as i32) as isize);
                y0r = x0r + x2r;
                y0i = x0i + x2i;
                y2r = x0r - x2r;
                y2i = x0i - x2i;
                y1r = x1r - x3i;
                y1i = x1i + x3r;
                y3r = x1r + x3i;
                y3i = x1i - x3r;
                x0r = *a.offset(j4 as isize) + *a.offset(j5 as isize);
                x0i =
                    *a.offset((j4 + 1 as i32) as isize) +
                        *a.offset((j5 + 1 as i32) as isize);
                x1r = *a.offset(j4 as isize) - *a.offset(j5 as isize);
                x1i =
                    *a.offset((j4 + 1 as i32) as isize) -
                        *a.offset((j5 + 1 as i32) as isize);
                x2r = *a.offset(j6 as isize) + *a.offset(j7 as isize);
                x2i =
                    *a.offset((j6 + 1 as i32) as isize) +
                        *a.offset((j7 + 1 as i32) as isize);
                x3r = *a.offset(j6 as isize) - *a.offset(j7 as isize);
                x3i =
                    *a.offset((j6 + 1 as i32) as isize) -
                        *a.offset((j7 + 1 as i32) as isize);
                y4r = x0r + x2r;
                y4i = x0i + x2i;
                y6r = x0r - x2r;
                y6i = x0i - x2i;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                x2r = x1r + x3i;
                x2i = x1i - x3r;
                y5r = wn4r * (x0r - x0i);
                y5i = wn4r * (x0r + x0i);
                y7r = wn4r * (x2r - x2i);
                y7i = wn4r * (x2r + x2i);
                x0r = y1r + y5r;
                x0i = y1i + y5i;
                *a.offset(j1 as isize) = wk1r * x0r - wk1i * x0i;
                *a.offset((j1 + 1 as i32) as isize) =
                    wk1r * x0i + wk1i * x0r;
                x0r = y1r - y5r;
                x0i = y1i - y5i;
                *a.offset(j5 as isize) = wk5r * x0r - wk5i * x0i;
                *a.offset((j5 + 1 as i32) as isize) =
                    wk5r * x0i + wk5i * x0r;
                x0r = y3r - y7i;
                x0i = y3i + y7r;
                *a.offset(j3 as isize) = wk3r * x0r - wk3i * x0i;
                *a.offset((j3 + 1 as i32) as isize) =
                    wk3r * x0i + wk3i * x0r;
                x0r = y3r + y7i;
                x0i = y3i - y7r;
                *a.offset(j7 as isize) = wk7r * x0r - wk7i * x0i;
                *a.offset((j7 + 1 as i32) as isize) =
                    wk7r * x0i + wk7i * x0r;
                *a.offset(j as isize) = y0r + y4r;
                *a.offset((j + 1 as i32) as isize) = y0i + y4i;
                x0r = y0r - y4r;
                x0i = y0i - y4i;
                *a.offset(j4 as isize) = wk4r * x0r - wk4i * x0i;
                *a.offset((j4 + 1 as i32) as isize) =
                    wk4r * x0i + wk4i * x0r;
                x0r = y2r - y6i;
                x0i = y2i + y6r;
                *a.offset(j2 as isize) = wk2r * x0r - wk2i * x0i;
                *a.offset((j2 + 1 as i32) as isize) =
                    wk2r * x0i + wk2i * x0r;
                x0r = y2r + y6i;
                x0i = y2i - y6r;
                *a.offset(j6 as isize) = wk6r * x0r - wk6i * x0i;
                *a.offset((j6 + 1 as i32) as isize) =
                    wk6r * x0i + wk6i * x0r;
                j += 2 as i32
            }
            k += m
        }
    };
}
unsafe extern "C" fn rftfsub(n: i32, a: *mut smpl_t,
                             nc: i32, c: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut kk: i32 = 0;
    let mut ks: i32 = 0;
    let mut m: i32 = 0;
    let mut wkr: smpl_t = 0.;
    let mut wki: smpl_t = 0.;
    let mut xr: smpl_t = 0.;
    let mut xi: smpl_t = 0.;
    let mut yr: smpl_t = 0.;
    let mut yi: smpl_t = 0.;
    m = n >> 1 as i32;
    ks = 2 as i32 * nc / m;
    kk = 0 as i32;
    j = 2 as i32;
    while j < m {
        k = n - j;
        kk += ks;
        wkr = 0.5f64 as smpl_t - *c.offset((nc - kk) as isize);
        wki = *c.offset(kk as isize);
        xr = *a.offset(j as isize) - *a.offset(k as isize);
        xi =
            *a.offset((j + 1 as i32) as isize) +
                *a.offset((k + 1 as i32) as isize);
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        let ref mut fresh16 = *a.offset(j as isize);
        *fresh16 -= yr;
        let ref mut fresh17 = *a.offset((j + 1 as i32) as isize);
        *fresh17 -= yi;
        let ref mut fresh18 = *a.offset(k as isize);
        *fresh18 += yr;
        let ref mut fresh19 = *a.offset((k + 1 as i32) as isize);
        *fresh19 -= yi;
        j += 2 as i32
    };
}
unsafe extern "C" fn rftbsub(n: i32, a: *mut smpl_t,
                             nc: i32, c: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut kk: i32 = 0;
    let mut ks: i32 = 0;
    let mut m: i32 = 0;
    let mut wkr: smpl_t = 0.;
    let mut wki: smpl_t = 0.;
    let mut xr: smpl_t = 0.;
    let mut xi: smpl_t = 0.;
    let mut yr: smpl_t = 0.;
    let mut yi: smpl_t = 0.;
    *a.offset(1 as i32 as isize) =
        -*a.offset(1 as i32 as isize);
    m = n >> 1 as i32;
    ks = 2 as i32 * nc / m;
    kk = 0 as i32;
    j = 2 as i32;
    while j < m {
        k = n - j;
        kk += ks;
        wkr = 0.5f64 as smpl_t - *c.offset((nc - kk) as isize);
        wki = *c.offset(kk as isize);
        xr = *a.offset(j as isize) - *a.offset(k as isize);
        xi =
            *a.offset((j + 1 as i32) as isize) +
                *a.offset((k + 1 as i32) as isize);
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        let ref mut fresh20 = *a.offset(j as isize);
        *fresh20 -= yr;
        *a.offset((j + 1 as i32) as isize) =
            yi - *a.offset((j + 1 as i32) as isize);
        let ref mut fresh21 = *a.offset(k as isize);
        *fresh21 += yr;
        *a.offset((k + 1 as i32) as isize) =
            yi - *a.offset((k + 1 as i32) as isize);
        j += 2 as i32
    }
    *a.offset((m + 1 as i32) as isize) =
        -*a.offset((m + 1 as i32) as isize);
}
unsafe extern "C" fn dctsub(n: i32, a: *mut smpl_t,
                            nc: i32, c: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut kk: i32 = 0;
    let mut ks: i32 = 0;
    let mut m: i32 = 0;
    let mut wkr: smpl_t = 0.;
    let mut wki: smpl_t = 0.;
    let mut xr: smpl_t = 0.;
    m = n >> 1 as i32;
    ks = nc / n;
    kk = 0 as i32;
    j = 1 as i32;
    while j < m {
        k = n - j;
        kk += ks;
        wkr = *c.offset(kk as isize) - *c.offset((nc - kk) as isize);
        wki = *c.offset(kk as isize) + *c.offset((nc - kk) as isize);
        xr = wki * *a.offset(j as isize) - wkr * *a.offset(k as isize);
        *a.offset(j as isize) =
            wkr * *a.offset(j as isize) + wki * *a.offset(k as isize);
        *a.offset(k as isize) = xr;
        j += 1
    }
    let ref mut fresh22 = *a.offset(m as isize);
    *fresh22 *= *c.offset(0 as i32 as isize);
}
unsafe extern "C" fn dstsub(n: i32, a: *mut smpl_t,
                            nc: i32, c: *mut smpl_t) {
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut kk: i32 = 0;
    let mut ks: i32 = 0;
    let mut m: i32 = 0;
    let mut wkr: smpl_t = 0.;
    let mut wki: smpl_t = 0.;
    let mut xr: smpl_t = 0.;
    m = n >> 1 as i32;
    ks = nc / n;
    kk = 0 as i32;
    j = 1 as i32;
    while j < m {
        k = n - j;
        kk += ks;
        wkr = *c.offset(kk as isize) - *c.offset((nc - kk) as isize);
        wki = *c.offset(kk as isize) + *c.offset((nc - kk) as isize);
        xr = wki * *a.offset(k as isize) - wkr * *a.offset(j as isize);
        *a.offset(k as isize) =
            wkr * *a.offset(k as isize) + wki * *a.offset(j as isize);
        *a.offset(j as isize) = xr;
        j += 1
    }
    let ref mut fresh23 = *a.offset(m as isize);
    *fresh23 *= *c.offset(0 as i32 as isize);
}
