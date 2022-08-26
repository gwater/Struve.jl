#              Struve Functions of first and second kind
#                     struveh(v, x), struvek(v, x)
#
# The Struve function is computed from its power series [1] for small arguments.
# Large arguments are computed from the asymptotic expansions for struvek(v, x) [2].
# Connection formula using bessely(v, x) [3] are used to compute struveh(v, x) from struvek(v, x) and conversely (struvek -> struveh).
# The power series form [1] can only be used for small x and when nu > x otherwise significant cancellation occurs.
# There is a region where the large argument expansion is not valid and where the power series is prone to cancellation.
# This region is also hard to indiscriminately use the recurrence relation because it is not always stable in either direction.
# For low orders (v < 30) it is stable to use forward recurrence for values of x less than 50.
# For the region of x ∈ (6, 39) we use Chebyshev approximation for v ∈ (0, 2) and then use forward recurrence to fill higher orders.
# Forward recurrence is also used starting with large argument asymptotic expansion up until order 30.
# When v > 30 and x > 30 (and where large argument expansion is not valid) we fall back to computing struvek(v, x)
# from its expansion in series of Bessel functions [4]. This method in general works best for nu > x but works well as long as x is not much larger than nu.
# Backward recurrence may work in this region however, the power series validity depends on x^2 so we would have to shift to very high values
# and use backward recurrence over a large number of orders which is not very efficient.
# If the arguments are not real we use the integral representation [5] using adaptive gaussian integration (QuadGK.jl).
# TODO: Investigate the large order expansions.
#
# [1] http://dlmf.nist.gov/11.2.E1
# [2] http://dlmf.nist.gov/11.6.E1
# [3] http://dlmf.nist.gov/11.2.E5
# [4] http://dlmf.nist.gov/11.4.E19
# [5] http://dlmf.nist.gov/11.5.E2
#

"""
struveh(nu, x)

Struve function of the first kind of order `nu`.

External links: [DLMF](http://dlmf.nist.gov/11.2.E1), [Wikipedia](https://en.wikipedia.org/wiki/Struve_function)
"""
struveh(v::Real, x::Real) = _struveh(float(v), float(x))
struveh(v, x) = _H_integral(v, x)

_struveh(v, x::AbstractFloat) = _H_integral(v, x)
_struveh(v::Float16, x::Float16) = Float16(_struveh(Float32(v), Float32(x)))

function _struveh(v, x::T) where T <: Union{Float32, Float64}
    if struvek_large_arg_cutoff(v, x)
        # use http://dlmf.nist.gov/11.2.E6
        return struvek_large_argument(v, x) + bessely(v, x)
    elseif struveh_power_series_cutoff(v, x)
        return struveh_power_series(v, x)
    elseif v < 30
        struvek_chebyshev_cutoff(x) && return struvek_chebyshev(v, x)[1] + bessely(v, x)
        v_floor = modf(v)[1]
        k0, k1 = struvek_large_argument(v_floor, x), struvek_large_argument(v_floor + 1, x)
        return struvek_up_recurrence(x, k1, k0, v_floor + 1, v)[1] + bessely(v, x)
    else
        return struveh_bessel_series(v, x)
    end
end

"""
struvek(nu, x)

Struve function of the second kind of order `nu`.

External links: [DLMF](http://dlmf.nist.gov/11.2.E5), [Wikipedia](https://en.wikipedia.org/wiki/Struve_function)
"""
struvek(v::Real, x::Real) = _struvek(float(v), float(x))
struvek(v, x) = _K_integral(v, x)

_struvek(v, x::AbstractFloat) = _K_integral(v, x)
_struvek(v::Float16, x::Float16) = Float16(_struvek(Float32(v), Float32(x)))

function _struvek(v::Real, x::T) where T <: Union{Float32, Float64}
    if struvek_large_arg_cutoff(v, x)
        return struvek_large_argument(v, x)
    elseif struveh_power_series_cutoff(v, x)
        return struveh_power_series(v, x) - bessely(v, x)
    elseif v < 30
        struvek_chebyshev_cutoff(x) && return struvek_chebyshev(v, x)[1]
        v_floor = modf(v)[1]
        k0, k1 = struvek_large_argument(v_floor, x), struvek_large_argument(v_floor + 1, x)
        return struvek_up_recurrence(x, k1, k0, v_floor + 1, v)[1]
    else
        return struveh_bessel_series(v, x) - bessely(v, x)
    end
end

# H_{nu}(x) Struve function power series
# http://dlmf.nist.gov/11.2.E1
# can use for x < 5 || nu > -0.75 + 0.41x + 0.023x^2
# struve H can be computed with forward recurrence only when x < nu
# backward recurrence may be stable always?
function struveh_power_series(v, x::T) where T
    MaxIter = 50000
    S = promote_type(T, Float64)
    v, x = S(v), S(x)
    out = zero(S)
    three_halves = S(3) / 2

    a = (x/2)^(v+1)
    (isinf(a) || gamma_max(v)) && return T(log_struveh_power_series(v, x))
    a /= gamma(v + three_halves) * gamma(three_halves)
    iszero(a) && return a

    t2 = (x/2)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= -inv((v + i + three_halves) * (i + three_halves)) * t2
    end
    return T(out)
end
function log_struveh_power_series(v, x::T) where T
    MaxIter = 10000
    out = zero(T)
    three_halves = T(3) / 2
    a = x / 2 
    t2 = a*a
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= -inv((v + i + three_halves) * (i + three_halves)) * t2
    end
    return exp(v*log(x/2) + log(out) - loggamma(v + three_halves) - loggamma(three_halves))
end
struveh_power_series_cutoff(nu, x::Float64) = x < 6 || nu > evalpoly(x, (-0.75, 0.41, 0.023))
struveh_power_series_cutoff(nu, x::Float32) = x < 26 || nu > evalpoly(x, (-10.0f0, 0.1f0, 0.012f0))

# K_{nu}(x) using large argument expansion
# http://dlmf.nist.gov/11.6.E1
function struvek_large_argument(v, x::T) where T
    MaxIter = 50000
    S = promote_type(T, Float64)
    v, x = S(v), S(x)
    out = zero(S)
    one_half = one(S) / 2

    a = (x/2)^(v-1)
    (isinf(a) || gamma_max(v)) && return T(log_struvek_large_argument(v, x))
    a *= gamma(one_half) / gamma(v + one_half)
    iszero(a) && return a

    t2 = (x/2)^(-2)
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= (i + one_half) * t2 * (v - one_half - i)
    end
    return T(out / π)
end

function log_struvek_large_argument(v, x::T) where T
    MaxIter = 50000
    out = zero(T)
    one_half = one(T) / 2
    a = 2 / x
    iszero(a) && return a
    t2 = a*a
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= (i + one_half) * t2 * (v - one_half - i)
    end
    return exp(log(out / π) + loggamma(one_half) - loggamma(v + one_half) + v*log(x/2)) 
end

struvek_large_arg_cutoff(nu, x) = x > 39 && nu < 1.23*x - 35
struvek_large_arg_cutoff(nu, x::Float32) = x > 17 && nu < 1.23*x - 23

#####
#####  Chebyshev approximation for struvek_{nu}(x)
#####

"""
    struvek_chebyshev(nu, x::T) where T <: Float64
Computes ``struvek_{nu}(x)`` for small to medium arguments x ∈ (6, 36) for any positive order using a Chebyshev approximation.
Forward recurrence is used to fill orders starting at low orders ν ∈ (0, 2). Two sets of Chebyshev coefficients are used for
x ∈ (6, 21) and for x ∈ (21, 36).
"""
function struvek_chebyshev(v, x)
    v_floor, _ = modf(v)
    k0, k1 = struvek_chebyshev_low_orders(v_floor, x)
    return struvek_up_recurrence(x, k1, k0, v_floor + 1, v)
end

# compute struvek for x ∈ (6, 39) and ν ∈ (0, 2) using chebyshev approximation in three regions
# optimized to return both (ν, ν + 1) in around the same time, therefore ν must be in (0, 1)
# no checks are performed on arguments
function struvek_chebyshev_low_orders(v, x)
    # need to rescale inputs according to
    #x0 = (x - lb) * 2 / (ub - lb) - 1
    v1 = v - 1
    v2 = v
    if x <= 12
        x1 = (x - 6) * 2 / 6 - 1
        P = clenshaw_chebyshev.(x1, struvek_cheb_weights_branch1)
        return clenshaw_chebyshev(v1, P), clenshaw_chebyshev(v2, P)
    elseif x <= 21
        x1 = (x - 12) * 2 / 9 - 1
        P2 = clenshaw_chebyshev.(x1, struvek_cheb_weights_branch2)
        return clenshaw_chebyshev(v1, P2), clenshaw_chebyshev(v2, P2)
    else
        x1 = (x - 21) * 2 / 18 - 1
        P3 = clenshaw_chebyshev.(x1, struvek_cheb_weights_branch3)
        return clenshaw_chebyshev(v1, P3), clenshaw_chebyshev(v2, P3)
    end
end

# only implemented for Float64 so far
struvek_chebyshev_cutoff(x) = (x <= 39.0 && x >= 6.0)

# uses the Clenshaw algorithm to recursively evaluate a linear combination of Chebyshev polynomials
function clenshaw_chebyshev(x, c)
    x2 = 2x
    c0 = c[end-1]
    c1 = c[end]
    for i in length(c)-2:-1:1
        c0, c1 = c[i] - c1, c0 + c1 * x2
    end
    return c0 + c1 * x
end

# Chebyshev coefficients for low orders and small arguments nu ∈ (0, 2), x ∈ (6, 12)
#=
using ArbNumerics, FastChebInterp
g(x) = struveh_power_series(BigFloat(x[1]), BigFloat(x[2])) - bessely(ArbFloat(BigFloat(x[1])), ArbFloat(BigFloat(x[2])))
lb, ub = [0,6], [2, 12]; # lower and upper bounds of the domain, respectively
x = chebpoints((15,18), lb, ub);
c = chebinterp(g.(x), lb, ub);
map(i->tuple(c.coefs[i,:]...), 1:size(c.coefs)[1])
=#
const struvek_cheb_weights_branch1 = (
    (0.8414037738679299, 0.14317176935040615, 0.0019505115785243532, -0.0005691582221003979, 0.00011414462507084262, -2.141715827764053e-5, 3.928056373983627e-6, -7.129126546948245e-7, 1.286264900503637e-7, -2.311782654689725e-8, 4.143186976718107e-9, -7.408746725469307e-10, 1.3223073993786884e-10, -2.356141390069198e-11, 4.192050321711585e-12, -7.448497580037278e-13, 1.3232707111135957e-13, -2.417490076435637e-14, 4.157990824077623e-15),
    (0.9547096033961995, 0.28861404449263034, -0.0032673414121821625, 0.00014781075728455736, 8.83567435808672e-6, -5.727180361953503e-6, 1.5927578130740753e-6, -3.6837911941122645e-7, 7.855852310661183e-8, -1.6002446439560773e-8, 3.1645533843759707e-9, -6.129024065678706e-10, 1.1687964834655522e-10, -2.202284490157776e-11, 4.110060834216697e-12, -7.610689727609875e-13, 1.4019557120153396e-13, -2.6447979608848174e-14, 4.710635174113257e-15),
    (0.19274534532054857, 0.14864106919743877, 0.0017012688968629608, -0.0002764872060537484, 4.544249851057666e-5, -7.726636721542653e-6, 1.3506857056745082e-6, -2.411153699621416e-7, 4.369322457957086e-8, -7.997643281788248e-9, 1.4728650399854992e-9, -2.721049892776596e-10, 5.0322689225544716e-11, -9.30271115657672e-12, 1.7173643442761686e-12, -3.164315774427127e-13, 5.82232307634035e-14, -1.1019991503686738e-14, 1.9605716227454615e-15),
    (0.0025874393078261547, 0.028195280261165356, 0.003906483781776456, -0.00035366913180259795, 3.778838440024789e-5, -4.388673531317863e-6, 5.423687993432057e-7, -7.134544490556748e-8, 1.0074025750294606e-8, -1.537403620871308e-9, 2.5306131712779807e-10, -4.437626248935845e-11, 8.143417216275482e-12, -1.536986022392129e-12, 2.943822551715542e-13, -5.670751933871964e-14, 1.0927267321444897e-14, -2.1743923534140103e-15, 3.945917666688577e-16),
    (-0.003042846141870137, -0.0003338360999496057, 0.0011457287377530164, -1.3062753141345805e-5, -5.562318921977263e-6, 1.31789702795212e-6, -2.318321127052029e-7, 3.697844179276834e-8, -5.606378064771359e-9, 8.193868396031999e-10, -1.1555567279556729e-10, 1.5590075924495412e-11, -1.967823435064828e-12, 2.1998985777603905e-13, -1.8126911889855663e-14, -1.753347659961947e-16, 5.780612614095752e-16, -2.0825034084013425e-16, 4.5950897408096754e-17),
    (-0.00010633338599974852, -0.0005759999096383057, 2.2540836623484335e-5, 2.895709521844765e-5, -3.735724321593328e-6, 3.8480851698695975e-7, -3.199478010023308e-8, 1.3007272397271832e-9, 2.811167271911154e-10, -1.0849517846633905e-10, 2.5813978002826575e-11, -5.295045590915662e-12, 1.0145778507454473e-12, -1.8752009653499358e-13, 3.396698141102735e-14, -6.1185977837827806e-15, 1.0665176303831853e-15, -1.9914625504213745e-16, 2.5127790794612196e-17),
    (2.8912643703330436e-5, -1.0208414039982949e-5, -2.3992888111514224e-5, 3.098624368056638e-6, 2.3345709468366125e-7, -9.395473444433735e-8, 1.8333649021822243e-8, -2.978113074160245e-9, 4.4401516928352756e-10, -6.276304894640289e-11, 8.514388687492209e-12, -1.111943777175078e-12, 1.3938317771057787e-13, -1.6621754104581202e-14, 1.8535956412416966e-15, -1.7556505856875782e-16, 4.144141948327742e-17, 1.8999036847367745e-17, 1.41275558638463e-17),
    (3.160146461539015e-7, 6.770385748705868e-6, -8.84767391499818e-7, -5.870192695164897e-7, 1.3621795000137387e-7, -1.3786509918873852e-8, 5.391485594185784e-10, 1.4921804027008386e-10, -5.330580799979856e-11, 1.1953069263548332e-11, -2.283837219591695e-12, 4.016326998861076e-13, -6.705821597596395e-14, 1.0782537748338806e-14, -1.6757284475630905e-15, 2.599118716354382e-16, -4.6478743243553127e-17, -3.518035271535372e-18, 8.954706832009907e-19),
    (-1.8324973754750047e-7, 1.527920058365505e-8, 3.072817171653234e-7, -5.7689358307697896e-8, -5.465658791102753e-9, 3.0251854538909887e-9, -6.077545710133883e-10, 9.16226778359545e-11, -1.1441149112745778e-11, 1.1333350410410767e-12, -5.888594708106496e-14, -1.0261363650414336e-14, 4.633703400404567e-15, -1.2013698007509953e-15, 2.7636940656890314e-16, -5.684453167387116e-17, 1.1680270626730411e-17, 7.266295274569112e-18, 1.2721305490496584e-18),
    (6.243968255457977e-9, -5.515275815718242e-8, 4.897411965563784e-9, 8.401174409777174e-9, -2.2571484881455044e-9, 2.0180695904341714e-10, 1.5611046650335596e-11, -9.670196237548094e-12, 2.38053600102002e-12, -4.581589223233378e-13, 7.770257566927765e-14, -1.210558117921416e-14, 1.7467123426663358e-15, -2.5063322770167584e-16, 3.233423983894695e-17, -9.386480308293255e-18, -8.121256312959003e-18, -7.700762394053705e-18, 4.979298866229724e-18),
    (4.5031176467261745e-10, 1.5440100694026544e-9, -2.7785920427438995e-9, 4.4923735421202444e-10, 1.2539207359888315e-10, -5.504494823035031e-11, 1.0409288287937051e-11, -1.273945078172784e-12, 7.886114913796261e-14, 1.1390216830419885e-14, -5.664581882067974e-15, 1.4585341486788097e-15, -2.9946788361858134e-16, 5.986898848095804e-17, -1.0779901093904798e-17, 3.8523514799733155e-18, -5.143170232620986e-18, 2.911285572028084e-18, 5.6298201697210265e-18),
    (-6.24968281801273e-11, 2.7725542548971617e-10, 4.5790645604887515e-11, -8.759880441745886e-11, 2.062243625435878e-11, -5.743238657506244e-13, -7.059848613592501e-13, 2.270741582933742e-13, -4.636704878325907e-14, 7.525693232492497e-15, -1.022640649524495e-15, 1.0098355864281444e-16, 9.603722299150536e-18, -8.049164342414597e-18, 1.3166571855322776e-17, -7.619333901246108e-18, -6.166120917688318e-18, -1.2488749261596412e-17, -2.6369202282108097e-18),
    (3.1215770567997312e-12, -1.9076246746191166e-11, 1.6240192219554057e-11, -5.690581915776475e-13, -1.829595954435418e-12, 5.97403697167055e-13, -9.048878057731838e-14, 4.150119389554837e-15, 1.8753871316450744e-15, -7.238332473911517e-16, 1.9699368113265635e-16, -4.0919510883148707e-17, 8.197271147620728e-18, -1.0446084930094262e-17, -4.690549795139961e-18, -1.794425232681166e-17, 1.0804690664772471e-17, 5.827719561760223e-18, -4.26006123060265e-18),
    (1.2480870015870587e-13, -5.056383831724078e-13, -8.61438505976208e-13, 5.963754729496639e-13, -8.929401312194873e-14, -1.9522344084558006e-14, 1.124444755494374e-14, -2.6879337850027665e-15, 4.3664431403764185e-16, -4.2038734569044894e-17, -2.7072897179807766e-18, 4.120488838006777e-18, 5.820659829896484e-18, 1.4435990687233315e-17, 6.5411179648065425e-18, 3.841737526351985e-18, 3.569875182552813e-18, 1.8291657053143335e-18, 2.923029134838929e-18),
    (-2.939253618015634e-14, 1.0407663361160601e-13, -4.7140004364555376e-14, -1.9432639796863764e-14, 1.5568632820767158e-14, -3.638457612387882e-15, 2.2053154055987376e-16, 8.838014903477148e-17, -6.494403193841361e-17, -9.537766666991001e-18, 5.449645376548901e-18, 1.0604041859680771e-17, 2.8009359235262556e-17, 1.2880747823202199e-17, -2.9224585439973445e-18, -4.007708117712367e-19, -1.7275772204233778e-18, 1.406689189475388e-19, 4.65966000352659e-18),
    (1.2782984547420206e-15, -2.90010071787261e-15, 5.791364714104922e-15, -2.2597648189193838e-15, -1.1591481006278994e-16, 2.764816107895135e-16, -9.49857476623745e-17, 1.2433872154253724e-17, -3.756235919016888e-18, 6.784696261598179e-18, 2.4096000721831817e-18, -4.564489374368292e-18, -8.635067969306773e-18, -3.879870215541777e-18, -3.961240756914403e-18, 1.2179984622678339e-17, 5.564790673110196e-18, 3.2039998188005063e-18, 5.2427198385076835e-18)
)

# Chebyshev coefficients for low orders and medium arguments nu ∈ (0, 2), x ∈ (12, 21)
#=
using ArbNumerics, FastChebInterp
g(x) = struveh_power_series(BigFloat(x[1]), BigFloat(x[2])) - bessely(ArbFloat(BigFloat(x[1])), ArbFloat(BigFloat(x[2])))
lb, ub = [0,12], [2, 21]; # lower and upper bounds of the domain, respectively
x = chebpoints((15,17), lb, ub);
c = chebinterp(g.(x), lb, ub);
map(i->tuple(c.coefs[i,:]...), 1:size(c.coefs)[1])
=#
const struvek_cheb_weights_branch2 = (
    (1.2053641958978805, 0.21858775256444044, -0.0008071349642947239, -6.330378806161045e-5, 1.6899870041925155e-5, -2.905434708373683e-6, 4.5161875725048586e-7, -6.753231583981462e-8, 9.912404414137635e-9, -1.4401721389243015e-9, 2.0795674652641251e-10, -2.990810478377471e-11, 4.2894604307523295e-12, -6.138297726863619e-13, 8.772022460264614e-14, -1.2462906523784808e-14, 1.7693907341477005e-15, -2.2291536808160006e-16),
    (1.6427710865468372, 0.39874507357639216, -0.0035748800401840418, 0.00022373913068511237, -1.7399951806623135e-5, 1.382133080109445e-6, -9.481019824813643e-8, 2.7563118577973624e-9, 8.312638624740248e-10, -2.6498522716171484e-10, 5.589018990383279e-11, -1.025678643224466e-11, 1.7557039746843088e-12, -2.8820693448885667e-13, 4.603371061246626e-14, -7.14700629832716e-15, 1.2260345244488003e-15, -1.1145768404080003e-16),
    (0.5759478567056728, 0.2338176432025742, 0.0006144440878328562, -9.170432783424626e-5, 1.1903503283365683e-5, -1.5550253220835535e-6, 2.06087353759185e-7, -2.7681837859191844e-8, 3.762842079845031e-9, -5.169548763826138e-10, 7.169828324483062e-11, -1.0028021017042137e-11, 1.412936199608981e-12, -2.0033746924457925e-13, 2.854376593474117e-14, -4.056232474086381e-15, 6.408816832346002e-16, -1.6718652606120004e-16),
    (0.10823044284275865, 0.07738767577512289, 0.0032554267613610524, -0.0002462592475363499, 2.246367630469418e-5, -2.2351388039352004e-6, 2.3352284942954804e-7, -2.5167574406584117e-8, 2.7725531130473444e-9, -3.1074339302602077e-10, 3.535916172009606e-11, -4.08355721291056e-12, 4.790856996163055e-13, -5.723003600161525e-14, 6.9680916810009245e-15, -8.923145446938268e-16, 1.0449157878825003e-16, -2.7864421010200007e-17),
    (0.009732666952917465, 0.013802775700677164, 0.0017437105599074439, -7.456927940937819e-5, 3.3333185187365222e-6, -5.8861522837885826e-8, -1.8436702767921826e-8, 4.39645230610059e-9, -7.334644495512274e-10, 1.0881538701766513e-10, -1.52664489942714e-11, 2.0748908831592707e-12, -2.7632541226701787e-13, 3.62861148739946e-14, -4.694835206872149e-15, 6.042007852641416e-16, -9.055936828315002e-17, -2.6122894697062506e-18),
    (5.2910538187603234e-5, 0.0010723509440241049, 0.00039260185328497064, 5.347353203297928e-6, -1.7665347253037992e-6, 2.341116222511015e-7, -2.694936491608646e-8, 2.9355622275822718e-9, -3.0819777242676203e-10, 3.115500460998431e-11, -2.9911618533734067e-12, 2.6371510371616e-13, -1.9505097709954847e-14, 7.627419682374067e-16, 9.114781369755e-17, -5.5239037744830094e-17, 8.272249987403127e-18, 1.4598888546847954e-17),
    (-6.142472263822157e-5, -2.9387028314054164e-5, 3.6217054824127404e-5, 5.7094545512276834e-6, -4.772572545519258e-7, 2.7316911665346842e-8, -4.0110839583440646e-10, -2.2003013365840026e-10, 5.175572200781783e-11, -8.653642348910814e-12, 1.2774818985740374e-12, -1.7691580133656382e-13, 2.3578561963227737e-14, -3.0974521907402977e-15, 3.707822617837535e-16, -5.224238797554467e-17, 2.755999404725898e-17, -1.4660114081294194e-17),
    (-3.1807717966322545e-6, -1.0449160016789362e-5, -6.770311548035437e-7, 8.727848061751831e-7, 1.0161836082247186e-8, -7.625308972003718e-9, 1.1740429356631713e-9, -1.4053360631357572e-10, 1.477638261259925e-11, -1.379087044556413e-12, 1.0728616528971591e-13, -4.943566049007699e-15, -4.390824205203028e-16, 1.733360143231223e-16, -5.714733986273522e-17, -4.686479834718188e-17, 3.704995188642622e-18, 1.3201330687648011e-17),
    (1.5933496660629197e-7, -2.810337685344888e-7, -3.6275681280589455e-7, 1.3907888204172754e-8, 1.3568179726449876e-8, -1.3266961506362183e-9, 6.240471864782963e-11, 4.348213025199429e-12, -1.740613239539354e-12, 3.252456200209006e-13, -4.949181759000151e-14, 6.800472051767182e-15, -8.762622580559434e-16, 1.2930115886568824e-16, 4.231407896023092e-18, -4.798205805649513e-18, -3.208069192920099e-18, -1.1554618917437674e-17),
    (1.4091842372737312e-8, 4.871079319831486e-8, -1.1487766876033792e-8, -7.566359747705563e-9, 9.909898510638615e-10, 8.003161288913848e-11, -2.4348972379612667e-11, 3.458721901837739e-12, -3.7416482210403563e-13, 3.1750229429608525e-14, -1.6164457948447773e-15, -1.2850652037858904e-16, 6.174435725572911e-17, -3.49078570491365e-17, -6.133899533994942e-18, -7.747079595011965e-18, 6.725454887990474e-18, 1.4802973661668754e-17),
    (-2.3081926386395386e-10, 2.160417386972096e-9, 1.8284203430433754e-9, -4.2273959475865174e-10, -9.02234344638548e-11, 2.3503510106641617e-11, -1.7417624751762875e-12, -3.7872476210019066e-14, 3.461204047013145e-14, -7.101190304896461e-15, 1.0987418223804632e-15, -1.3782378260856752e-16, 2.504688809222274e-17, 1.3687164587105353e-17, 3.046136921846746e-17, 1.9384801164032094e-17, 1.3575687694205004e-18, 9.675335151796132e-18),
    (-2.8250107558020578e-11, -1.778059784407929e-10, 9.061439774578732e-11, 3.69282227884566e-11, -1.1146032359349214e-11, -5.0260327376423986e-14, 2.777826126684704e-13, -4.978810781275265e-14, 5.6092895947454545e-15, -3.9171669808840734e-16, 1.1581730272422424e-17, 3.901345757532067e-18, 1.245808397816409e-17, -2.3467992502329146e-17, -1.3576480402471009e-17, 9.491229154135667e-18, -3.561680898975727e-19, 2.7375639817562456e-18),
    (1.9811480727741816e-14, -7.895707367607144e-12, -7.258773074353808e-12, 2.8467316070868392e-12, 3.713839974165249e-13, -1.9676542685840308e-13, 2.179151405041764e-14, 1.9513109601126715e-16, -4.56368541112302e-16, 1.0468887579698434e-16, -6.601953009075233e-18, -2.9700873090014583e-18, 9.607791561906403e-18, -4.891729463918932e-18, 2.6561543158880066e-17, 1.901222314724583e-18, 2.6164475133934568e-17, 3.207224152991544e-17),
    (-5.803066069795924e-15, 5.881296700293187e-13, -3.680740552688895e-13, -1.5036352000997716e-13, 6.64423860379044e-14, -2.0523974689191387e-15, -2.069636457527078e-15, 4.676091224695159e-16, -5.531645245646007e-17, -2.2700482561950138e-17, 2.0063398295062025e-17, 8.771486188555819e-18, -1.799163405634218e-17, -6.471858244377519e-18, -6.979536530335166e-18, 3.141509457548149e-17, 3.3274691274143066e-17, -6.118082439341655e-18),
    (2.1646848937533996e-15, 1.6794252638594138e-14, 2.599035398428528e-14, -1.1699846789626251e-14, -1.4863685779552956e-15, 1.135842670018254e-15, -1.3344986659075796e-16, -5.3034213897777136e-18, -1.524338815806058e-17, -1.5726088958709469e-19, -1.0823333824220964e-17, 6.971976519995742e-18, 1.3888980495313444e-17, -1.6499444717442752e-17, 5.473001236211586e-18, -1.268515036407275e-17, -2.2802480861532287e-19, -7.039235752027455e-18),
    (1.750233944703188e-16, -1.7189767483698836e-15, 9.46710533460137e-16, 5.828888453002618e-16, -2.711865915816354e-16, -1.0309320780426231e-17, 1.5757743923308922e-17, -2.8378037038838645e-18, -9.020434706059913e-18, -4.2335829154333153e-19, 5.0389532346599395e-18, 1.7066770041933533e-17, 9.329746612435746e-18, 1.8586640173942127e-17, -2.0566661082738725e-18, 1.858492466611196e-17, 1.8581740395233762e-17, -1.4802973661668754e-17)
)

# Chebyshev coefficients for low orders and medium arguments nu ∈ (0, 2), x ∈ (21, 38)
#=
using ArbNumerics, FastChebInterp
g(x) = struveh_power_series(BigFloat(x[1]), BigFloat(x[2])) - bessely(ArbFloat(BigFloat(x[1])), ArbFloat(BigFloat(x[2])))
lb, ub = [0,21], [2, 39]; # lower and upper bounds of the domain, respectively
x = chebpoints((15,17), lb, ub);
c = chebinterp(g.(x), lb, ub);
map(i->tuple(c.coefs[i,:]...), 1:size(c.coefs)[1])
=#
const struvek_cheb_weights_branch3 = (
    (1.8426809188379392, 0.4163151016442681, -0.0031455661257598095, 0.0001241577253354793, -1.5244679130608855e-6, -1.0063088658561698e-6, 2.6077493746323394e-7, -5.034008919207968e-8, 8.820748755163872e-9, -1.478238120804951e-9, 2.4183025816990793e-10, -3.89954984913716e-11, 6.230290197086073e-12, -9.892027790404144e-13, 1.563016328565583e-13, -2.474273509390104e-14, 3.873154520417801e-15, -4.458307361632001e-16),
    (2.7872501673820294, 0.7422571893475773, -0.006359403582013682, 0.0004392955575590576, -4.0584781489437266e-5, 4.2488721430342e-6, -4.736982334567531e-7, 5.445522294242181e-8, -6.317523280650738e-9, 7.258302670219042e-10, -8.074397858507819e-11, 8.379969383648545e-12, -7.442351703254416e-13, 3.894521551655823e-14, 4.6986924155422616e-15, -2.0941853915478444e-15, 3.343730521224001e-16, 0.0),
    (1.2828274995862812, 0.4727789998799774, 0.00016137615617208757, -7.521463228481614e-5, 1.1682602041296445e-5, -1.7059710758051966e-6, 2.4821991279483193e-7, -3.630192447751896e-8, 5.341907118622628e-9, -7.905991240878264e-10, 1.1761228609363194e-10, -1.757781199832181e-11, 2.638286789556632e-12, -3.9752655994055027e-13, 6.011248824419341e-14, -9.15834673893208e-15, 1.4210854715202005e-15, -2.2291536808160006e-16),
    (0.38452665500251243, 0.20128068481870462, 0.004945467123585029, -0.00040936737844719584, 4.161355786810126e-5, -4.672970495434197e-6, 5.55771574265904e-7, -6.85711147865239e-8, 8.675406881630438e-9, -1.1175681024143366e-9, 1.459187570312056e-10, -1.9251533347894683e-11, 2.5610515418963486e-12, -3.4298581812388173e-13, 4.623853043370042e-14, -6.2590455694161765e-15, 8.916614723264002e-16, -8.359326303060002e-17),
    (0.07871886868646577, 0.05747122075260976, 0.003724913963066685, -0.00021758281766563373, 1.6255757264094045e-5, -1.325793988934931e-6, 1.0880556404924578e-7, -8.21553836336375e-9, 4.5033715782894923e-10, 1.1579560662548141e-11, -1.0013093158356589e-11, 2.387071329207923e-12, -4.54008476353522e-13, 7.859178313587247e-14, -1.2900764334795677e-14, 2.039109621894871e-15, -3.622374731326001e-16, 2.0898315757650005e-17),
    (0.011304448418992783, 0.011204387562353882, 0.0013842626998823702, -3.744314628648039e-5, 1.5616397689921421e-7, 2.176655993802925e-7, -4.3361349384426164e-8, 6.745283582901641e-9, -9.677659531096168e-10, 1.338450365363049e-10, -1.8140097433091208e-11, 2.4263206276857513e-12, -3.2119686434886786e-13, 4.212361705829606e-14, -5.485858907661832e-15, 6.797394890964806e-16, -8.098097356089377e-17, -9.578394722256252e-18),
    (0.0011446134035160045, 0.0014995549872893206, 0.00031185928662801933, 4.870140462663733e-6, -1.2186357158134281e-6, 1.4724648345583144e-7, -1.571604620210656e-8, 1.567561443180927e-9, -1.4381350599648232e-10, 1.1191459734977566e-11, -5.027538549821011e-13, -5.774190571797999e-14, 2.3712136142559044e-14, -5.200193531844012e-15, 9.097327840664596e-16, -1.5148557789431977e-16, 7.728023014547659e-18, 7.074950647121095e-18),
    (7.936114463101824e-5, 0.00013429072617973995, 4.526664532979743e-5, 3.462489757843435e-6, -2.889601390205165e-7, 1.621062269668272e-8, -1.8485987769204041e-10, -1.4587407340901602e-10, 3.423214939414424e-11, -5.837114987330688e-12, 8.835495031855377e-13, -1.2548074983944218e-13, 1.715188039552987e-14, -2.1996163464830884e-15, 2.348216086447811e-16, -6.100784365709807e-17, 3.2694435394292295e-17, 1.056480611055679e-17),
    (3.305240059051955e-6, 7.210563384863732e-6, 4.185453280626156e-6, 7.267333530441801e-7, -1.4841861306921621e-8, -2.5935857576537622e-9, 5.309032773488351e-10, -7.093231722376564e-11, 8.029982918490384e-12, -7.936026269994414e-13, 6.377347926734285e-14, -2.6003658794481776e-15, -4.095226825811481e-16, 1.4172475718795867e-16, -3.327974512593851e-17, 4.926954813632167e-17, 4.8185345963836166e-18, -1.0394735181539456e-17),
    (2.3515920128816375e-8, 1.1616374313448238e-7, 2.1428813167369913e-7, 8.235264070870957e-8, 4.6588935999759974e-9, -7.837305489325846e-10, 6.570236502899559e-11, -2.7500483460733824e-12, -3.126032881711072e-13, 1.1065067046619782e-13, -2.123415346937319e-14, 3.408123408457798e-15, -5.008398726614593e-16, 8.597082596861855e-17, -2.226197236877278e-17, -1.511454891834504e-17, 2.5211706477873817e-17, 2.0898315757650005e-17),
    (-6.821661364653666e-9, -1.2829003109731098e-8, 1.0869170998849215e-9, 4.889655889588383e-9, 1.0165751532050328e-9, -5.06451530510389e-11, -3.6161989124621897e-12, 1.1897834031035531e-12, -1.847573077281464e-13, 2.229482999895978e-14, -2.1483018770017163e-15, 1.7959337163482078e-16, -2.3780443506473404e-20, -2.5139974309380734e-17, 3.2049158221253845e-17, 1.0386450763342468e-17, -2.3628151809399244e-17, 3.896665125645157e-17),
    (-4.4170708712046257e-10, -9.18962035841212e-10, -6.026260770889342e-10, 4.351995146683136e-11, 8.611373127591438e-11, 5.252928403513797e-12, -1.3103379560419413e-12, 1.3094954572625068e-13, -6.014591066735279e-15, -8.348692313216403e-16, 2.37726512384914e-16, -2.824534003001389e-17, 2.2817642634047198e-17, 4.9584292805257635e-17, 1.676260448400854e-17, 4.9255462736171145e-17, -1.0228389640444183e-17, -4.398034224388258e-18),
    (-7.1236628179751525e-12, -1.1429222893214948e-11, -3.17788287235429e-11, -1.4398515740988387e-11, 2.1060627752625567e-12, 1.0643879400778358e-12, -6.671504726954648e-14, -5.775964460855574e-15, 1.950600093847615e-15, -3.1634846252453553e-16, 1.0984420785998259e-16, 3.2643896339189964e-17, 2.4464710182849977e-17, 3.5472920631676756e-17, 1.2570407436968895e-16, -2.001080097723498e-17, -5.67290753996727e-17, 6.946121918390465e-18),
    (6.068474528714394e-13, 1.2484033006826216e-12, 2.0302353851679809e-13, -8.199321152288225e-13, -2.0752139989773238e-13, 5.5695594843059837e-14, 6.641692439463549e-15, -1.6719502353545431e-15, 1.3389862632129946e-16, -6.099037971169688e-17, -2.429135004616317e-17, -2.353592754101515e-17, -3.779240215796303e-17, -2.0132087858904546e-17, 1.6257174904733667e-17, 1.4987861445465922e-17, -2.222254810781979e-17, -2.9626355834819587e-18),
    (3.5988695683683736e-14, 5.3944161580493486e-14, 7.157977619599883e-14, 6.53687440194854e-15, -1.6667247722629775e-14, -1.3727465863347025e-15, 8.302766576007026e-16, -8.640323759459079e-17, -9.818538374862911e-18, 6.0588136910083824e-18, 2.349719321040363e-17, -9.706543256395945e-18, -8.162466881245411e-18, 3.23693347036566e-17, 1.9783474070421363e-17, -2.5709504849124885e-17, 2.2073227371358687e-18, 1.0285889786968363e-17),
    (-1.4236977609899067e-16, -3.897300723909553e-16, 1.3345083760394254e-15, 2.0432192036508327e-15, -3.769778699930298e-17, -2.450996938363859e-16, 5.117235947138352e-18, -2.1366994993435147e-17, 8.902356067671317e-18, -2.2698377644076436e-18, 1.56321661735351e-17, 2.641810633159867e-17, 9.830805457249974e-19, -8.226323443063988e-18, -9.72629231992767e-19, -1.0157054784950113e-18, 7.536357321126083e-18, 1.1755302613678128e-17)
)

function struvek_up_recurrence(x::T, knu, knum1, nu_start, nu_end) where T
    x_2 = x / 2
    two_x = 2 / x
    c = x_2^nu_start / (sqrt(π) * gamma(nu_start + T(3)/2)) 
    while nu_start < nu_end + 0.5 # avoid inexact floating points when nu isa float
        knum1, knu = knu, muladd(nu_start*two_x, knu, c - knum1)
        c *= x_2 / (nu_start + T(3)/2)
        nu_start += 1
    end
    return knum1, knu
end

# Expansion of struveh in series of Bessel functions
# http://dlmf.nist.gov/11.4.E19
# The form given by 11.4.E20 is more prone to cancellation
# In general, this is accurate when nu > x/2 though more terms are needed when x > nu
# A simple translation of the formula given by 11.4.E19 will be very slow as it requires calls to besselj
# each time within the loop
# Because we only require besselj(k, x) and then besselj(k+1, x) during the loop we can use recurrence
# Unfortunately, forward recurrence is only stable when nu < x but we usually employ this when nu > x so we can't use forward
# Backward recurrence for besselj is always stable so we will employ that
# The difficult part is we can no longer check for convergence as we sum through the loop as we need to know how many terms
# This was solved by defining regions of convergence for a set amount of terms (e.g., 45, 75, 150)
# We need to be careful not to indiscriminately use a large amount of terms as besselj will underflow which will make recurrence useless
# These optimizations speed up the code by roughly 20x even when using more iterations within the sum 
function struveh_bessel_series(v, x::T) where T
    x2 = x / 2
    two_x = 2 / x
    out = zero(T)

    # need to be careful not to start loop too high as besselj -> 0 and could underflow
    if v > evalpoly(x, (-5.0, 0.001, 0.021))
        Iter = 45
    elseif v > evalpoly(x, (-3.0, 0.15, 0.008))
        Iter = 75
    else
        Iter = 150
    end

    # compute besselj(v, x) and besselj(v+1, x)
    jnup1 = besselj(Iter + T(3)/2 + v, x)
    jnu = besselj(Iter + T(1)/2 + v, x)

    # avoid overflow
    x2_pow = x2^(Iter/2)
    a = x2_pow / gamma(Iter + 1)
    for k in Iter:-1:0
        out += a / (k + T(1)/2) * jnu
        a *= k / x2
        jnup1, jnu = jnu, muladd((k + T(1)/2 + v)*two_x, jnu, -jnup1)
    end
    return out*sqrt(x2 / π) * x2_pow
end