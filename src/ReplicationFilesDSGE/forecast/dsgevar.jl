# This script holds functions mapping the solution to a DSGE model
# to the corresponding VAR system

function DSGE_to_VAR!(m::AbstractDSGEModel, observables::Vector{Symbol},
                      exogenous_shocks::Vector{Symbol}, lags::Int;
                      regime_switching::Bool = false, n_regimes::Int = 2,
                      regime::Int = 1, zero_DD::Bool = false,
                      MM::Matrix{S} =
                      zeros(length(observables), length(exogenous_shocks))) where {S<:Real}
    para = map(x -> x.value, m.parameters)
    return DSGE_to_VAR!(m, para, observables, exogenous_shocks, lags;
                        regime_switching = regime_switching,
                        n_regimes = n_regimes, regime = regime, zero_DD = zero_DD, MM = MM)
end

function DSGE_to_VAR!(m::AbstractDSGEModel, para::Vector{S},
                      observables::Vector{Symbol},
                      exogenous_shocks::Vector{Symbol}, lags::Int;
                      regime_switching::Bool = false, n_regimes::Int = 2,
                      regime::Int = 1,
                      zero_DD::Bool = false,
                      MM::Matrix{S} =
                      zeros(length(observables), length(exogenous_shocks)),
                      get_VAR::Bool = true) where {S<:Real}
    DSGE.update!(m, para)
    if regime_switching
        regime_system = compute_system(m; regime_switching = true, n_regimes = n_regimes)
        system = System(regime_system, regime)
        update_system!(m, system, observables, exogenous_shocks, zero_DD = zero_DD)
    else
        system = compute_system(m)
        update_system!(m, system, observables, exogenous_shocks, zero_DD = zero_DD)
    end

    return DSGE_to_VAR(system[:TTT], system[:RRR], system[:QQ], system[:DD],
                       system[:ZZ], system[:EE], MM, lags; get_VAR = get_VAR)
end

function DSGE_to_VAR(TTT::Matrix{S}, RRR::Matrix{S}, QQ::Matrix{S},
                     DD::Vector{S}, ZZ::Matrix{S}, EE::Matrix{S},
                     MM::Matrix{S}, nlags::Int; get_VAR::Bool = true) where {S<:Real}
## description:
## nlags: number VAR lags desired
## solution to DSGE model - delivers transition equation for the state variables  S_t
## transition equation: s_t = C + TTT s_{t-1} +  RRR ϵ_t, where var(ϵ_t) = QQ
## define the measurement equation: X_t = ZZ s_t + D + u_t
## where u_t = η_t + MM * ϵ_t with var(η_t) = EE
## where var(u_t) = HH = EE + MM QQ MM', cov(eps_t,u_t) = VV = QQ * MM'
## this is for no coint no constant


    nobs = size(ZZ,1)

    yyyyd = zeros(nobs,nobs)
    xxyyd = zeros(nlags*nobs,nobs)
    xxxxd = zeros(nlags*nobs,nlags*nobs)

    HH = EE+MM*QQ*MM';
    VV = QQ*MM';

    ## Compute nlags autocovariances

    ## Initialize Autocovariances
    GAMM0 = zeros((nobs)^2,nlags+1)

    GA0 =  DSGE.solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
    Gl   = ZZ*GA0*ZZ' + ZZ*RRR*VV + (ZZ*RRR*VV)' + HH
    GAMM0[:,1] = vec(Gl)

    TTl = copy(TTT)
    for l = 1:nlags
        Gl = ZZ*(TTl*GA0)*ZZ' + ZZ*(TTl*RRR*VV)
        GAMM0[:,l+1] = vec(Gl)
        TTl = TTl*TTT
    end

    ## Create limit cross product matrices

    yyyyd = reshape(GAMM0[:,1],nobs,nobs) + DD*DD'

    ## cointadd are treated as the first set of variables in XX
    ## coint    are treated as the second set of variables in XX
    ## composition: cointadd - coint - constant - lags
    yyxxd = zeros(nobs,nlags*nobs)
    xxxxd = zeros(nlags*nobs,nlags*nobs)
    for rr = 1:nlags;
        ## E[yy,x(lag rr)]
        yyxxd[:,nobs*(rr-1)+1:nobs*rr] =  reshape(GAMM0[:,rr+1],nobs,nobs) + DD*DD'


        ## E[x(lag rr),x(lag ll)]
        for ll = rr:nlags;
            yyyydrrll = reshape(GAMM0[:,ll-rr+1],nobs,nobs)+DD*DD';
            xxxxd[nobs*(rr-1)+1:nobs*rr,nobs*(ll-1)+1:nobs*ll] =  yyyydrrll
            xxxxd[nobs*(ll-1)+1:nobs*ll,nobs*(rr-1)+1:nobs*rr] =  yyyydrrll'
        end
    end

    xxyyd = convert(Matrix{S}, yyxxd')

    if get_VAR
        β = \(xxxxd, xxyyd)
        Σ = yyyyd - xxyyd' * β
        Σ += Σ'
        Σ ./= 2.
        return β, Σ
    else
        return yyyyd, xxyyd, xxxxd
    end
end
