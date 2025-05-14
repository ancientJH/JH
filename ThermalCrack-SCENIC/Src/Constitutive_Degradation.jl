function Degradation(ϕ::Float64, mat::Material)
    kmin = 1e-10;
    flagD = mat.flagD
    if flagD == QuadraticDegradation
        gd = (1.0 - ϕ)^2 + kmin; 
        gd′ = -2.0(1.0 - ϕ);
        gd′′ = 2.0;
    elseif flagD == WuDegradation
        a₁ = mat.a₁;
        fact = (1.0 - ϕ)^2 + a₁*ϕ*(1 - 0.5ϕ);
        gd = (1.0 - ϕ)^2/fact + kmin;
        gd′ = -a₁*(1.0 - ϕ)/fact^2;
        gd′′ = 2a₁^2/fact^3 - 3a₁/fact^2;
        gd′′ = gd′′>0 ? gd′′ : 0;
    elseif flagD == HughesDegradation
        s = mat.s;
        fact = ((s-1)/s)^((1.0 - ϕ)^2);
        fact2 = log(s/(s-1));
        gd = s*(1 - fact) + kmin;
        gd′ = -2s*(1.0 - ϕ)*fact*fact2;
        gd′′ = 2s*fact2*fact*(1-2fact2*(1.0 - ϕ)^2);
        gd′′ = gd′′>0 ? gd′′ : 0;
    end
    return gd, gd′, gd′′
end

function Constitutive(ε::SymmetricTensor{2,3,Float64}, mat::Material, gdn::Float64)
    Heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x, 0.5)))
    flag = mat.flag
    if flag == Isotropic
        μ = mat.G;
        λ = mat.λ;
        I = one(SymmetricTensor{2,3});
        𝕀 = one(SymmetricTensor{4,3});
        D = λ * I ⊗ I + 2μ * 𝕀;
        Ψ⁺ = 0.5 * ε ⊡ D ⊡ ε;
        Dᵉ = gdn * D; 
        σ = Dᵉ ⊡ ε;
    elseif flag == VolDev 
        K = mat.K;
        G = mat.G;
        I = one(SymmetricTensor{2,3});
        𝕀 = one(SymmetricTensor{4,3});
        D⁺ = K * Heaviside(tr(ε)) * I ⊗ I + 2G * (𝕀 - 1 / 3 * I ⊗ I);
        D⁻ = K * Heaviside(-tr(ε)) * I ⊗ I;
        Ψ⁺ = tr(ε) >0 ? 0.5*K*tr(ε)^2 + G*dev(ε) ⊡ dev(ε) : G*dev(ε) ⊡ dev(ε);
        σ⁺ = tr(ε) >0 ? K*tr(ε)*I + 2G*dev(ε) : 2G*dev(ε); 
        σ⁻ = tr(ε) <0 ? K*tr(ε)*I  : zero(Tensor{2, 3}); 
        Dᵉ = gdn * D⁺ + D⁻;
        σ = gdn *σ⁺ + σ⁻;
    elseif flag == Spectral 
        εₙ, Vₙ = eigen(ε);        # εₙ = eigvals(ε)
        μ = mat.G;
        λ = mat.λ;
        bracket₊(a::AbstractFloat) = a > 0 ? a : 0;
        bracket₋(a::AbstractFloat) = a < 0 ? a : 0;
        H₁₂⁺(x::AbstractFloat, y::AbstractFloat) = x ≠ y ? (bracket₊(x) - bracket₊(y)) / (x - y) : Heaviside(x);
        H₁₂⁻(x::AbstractFloat, y::AbstractFloat) = x ≠ y ? (bracket₋(x) - bracket₋(y)) / (x - y) : Heaviside(-x);
        I = one(SymmetricTensor{2,3});
        Ψ⁺ = tr(ε) >0 ? λ/2*tr(ε)^2 : 0.0;
        for e in εₙ;Ψ⁺ += e>0 ? μ*e^2 : 0.0;end;
        σ⁺ = λ*bracket₊(tr(ε))*I + 2μ*(bracket₊(εₙ[1])*Vₙ[:,1]⊗Vₙ[:,1] + bracket₊(εₙ[2])*Vₙ[:,2]⊗Vₙ[:,2] + bracket₊(εₙ[3])*Vₙ[:,3]⊗Vₙ[:,3]);
        σ⁻ = λ*bracket₋(tr(ε))*I + 2μ*(bracket₋(εₙ[1])*Vₙ[:,1]⊗Vₙ[:,1] + bracket₋(εₙ[2])*Vₙ[:,2]⊗Vₙ[:,2] + bracket₋(εₙ[3])*Vₙ[:,3]⊗Vₙ[:,3]);
        D⁺ = λ * Heaviside(tr(ε)) * I ⊗ I;
        D⁻ = λ * Heaviside(-tr(ε)) * I ⊗ I;
        for a in 1:3
            D⁺ += 2μ * Heaviside(εₙ[a]) * (Vₙ[:, a] ⊗ Vₙ[:, a]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, a]);
            D⁻ += 2μ * Heaviside(-εₙ[a]) * (Vₙ[:, a] ⊗ Vₙ[:, a]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, a]);
            for b in 1:3
                if a ≠ b
                    D⁺ += μ * H₁₂⁺(εₙ[a], εₙ[b]) * ((Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, b]) + (Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, b] ⊗ Vₙ[:, a]));
                    D⁻ += μ * H₁₂⁻(εₙ[a], εₙ[b]) * ((Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, a] ⊗ Vₙ[:, b]) + (Vₙ[:, a] ⊗ Vₙ[:, b]) ⊗ (Vₙ[:, b] ⊗ Vₙ[:, a]));
                end
            end
        end
        Dᵉ = gdn * D⁺ + D⁻;
        σ = gdn *σ⁺ + σ⁻;
    elseif flag == HybridSpectral 
        εₚ= eigvals(ε)
        μ = mat.G;
        λ = mat.λ;
        I = one(SymmetricTensor{2,3});
        𝕀 = one(SymmetricTensor{4,3});
        D = λ * I ⊗ I + 2μ * 𝕀;
        Dᵉ = gdn * D; 
        σ = Dᵉ ⊡ ε;
        Ψ⁺ = sum(εₚ)>0 ? λ/2*sum(εₚ)^2 : 0.0
        for e in εₚ
            Ψ⁺ += e>0 ? μ*e^2 : 0.0
        end
    end
    return Ψ⁺, Dᵉ, σ, σ⁺
end