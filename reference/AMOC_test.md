# Conduct an AMOC hypothesis test

Conduct an at-most one changepoint hypothesis test for changes in the
covariance operator of functional data based on the FKWC (functional
Kruskal–Wallis covariance changepoint) procedures outlined by Ramsay and
Chenouri (2025).

## Usage

``` r
amoc_test(data, ranks = NULL, depth = c("RPD", "FM", "LTR", "FMd", "RPDd"))
```

## Arguments

- data:

  Data in `matrix` or `data.frame` form, where each row is an
  observation and each column is a dimension.

- ranks:

  Optional if data is already ranked.

- depth:

  Depth function of choice.

## Value

A list consisting of:

- `$changepoint` : Index of the estimated changepoint.

- `$pvalue` : The p-value based on the null distribution.

- `$ranks` : A `vector` of depth-based ranks for each observation.

- `$method` : A `string` `"AMOC test (KWCChangepoint)"`

## Note

The options for the `depth` argument are as follows:

- `RPD`: Random projection depth

- `FM`: Frainman-Muniz depth

- `LTR`: \\L^2\\-root depth, most suitable for detecting changes in the
  norm

- `FMd`: Frainman-Muniz depth of the data and its first order derivative

- `RPDd`: Random projection depth of the data and its first order
  derivative

  The depth arguments that incorporate the first order derivative (which
  is approximated using
  [fda.usc::fdata.deriv](https://moviedo5.github.io/fda.usc/reference/fdata.deriv.html))
  result in a more robust detection of changes in the covariance
  structure (Ramsay and Chenouri, 2025).

## References

Ramsay, K., & Chenouri, S. (2025). Robust changepoint detection in the
variability of multivariate functional data. Journal of Nonparametric
Statistics. https://doi.org/10.1080/10485252.2025.2503891

## Examples

``` r
set.seed(11)
test_data <- rbind(replicate(3,rnorm(200,1,1)), #before changepoint
                   replicate(3,rnorm(200,1,5))) #after changepoint

amoc_test(test_data)
#> $changepoint
#> [1] 200
#> 
#> $pvalue
#> [1] 0
#> 
#> $ranks
#>   [1] 261 259 232 186 308 349 350 318 237 202 294 337 220 324 249 364 360 299
#>  [19] 339 357 272 372 387 380 356 321 395 325 289 334 309 330 238 258 231 394
#>  [37] 280 263 329 348 270 317 239 298 254 286 248 312 234 375 274 260 296 226
#>  [55] 326 389 365 366 221 253 307 195 400 228 320 290 215 303 384 292 371 247
#>  [73] 191 362 345 379 310 262 283 391 397 276 315 359 236 227 206 268 354 205
#>  [91] 331 347 374 244 243 266 314 291 222 201 265 279 370 218 198 346 355 332
#> [109] 264 242 287 328 376 233 323 378 295 383 340 393 194 327 187 277 368 240
#> [127] 358 197 319 363 361 316 219 342 382 386 257 208 217 229 216 336 304 388
#> [145] 353 390 306 302 293 275 301 385 300 190 284 373 267 196 224 333 311 281
#> [163] 273 367 245 399 271 255 396 250 377 297 278 209 210 351 200 241 252 230
#> [181] 344 343 235 341 282 352 269 313 392 305 223 213 381 214 369 251 288 225
#> [199] 398 338 199  53 145  63  61 147  44  69 175 166 180 161 207  25   4 116
#> [217] 162  23  55  38  36  96  58 137  79 104 109 179  81  21 131 193 164 119
#> [235]  50   7 103  78  16 121 155 335 154 125  20  48  89  97  60  54  35   8
#> [253]  26 105  86  92 152 120  34 163  11 138  72 141  87 113  77   1 112  68
#> [271]  18 177 181 159   2  90 178 107 140  99  93  94  39  73  31 174 142  85
#> [289]  46  13 150 184  51 108 158 160 183  84  70 149 167 171 246 144 101  80
#> [307]  28  88 139  91 127  42 156 129  66 165  43  59 111  74  40  75  12 151
#> [325] 168  95  83 146   3  64  14 134 153 172 122  57 212 143 188 157 133 169
#> [343]  45 211  15 204 170  47  30 118  19  71 136 106 135  17   9   5 110 123
#> [361]  82 203  32 256  29  24 114 115 322  33 185  37 126  27 117 176 102  62
#> [379]  65 189 148  98 124   6  67  76  10 132 128  56 285  49  22  52 100 182
#> [397] 192 173  41 130
#> 
#> $method
#> [1] "AMOC test (KWCChangepoint)"
#> 
```
