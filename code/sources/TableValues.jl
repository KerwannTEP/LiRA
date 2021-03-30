# Tables shifted by 1

# tables which don't depend on x,q

"""
    tabIln

Table of values of `I_ln=I(m,l,n)` used for computing Legendre integrals.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
const tabIln = zeros(Float64,2*N+1,2*N+1)

"""
    tabIaln

Table of values of `Ia_ln=Ia(m-1;l,m-1)` used for computing Legendre integrals.
See `Aoki et al. (1979)` for a definition (a=7/4).
"""
const tabIaln = zeros(Float64,2*N+1)

"""
    tabJln

Table of values of `J_ln=J(m,l,n)` used for computing Legendre integrals.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
const tabJln = zeros(Float64,2*N+1,2*N+1)

"""
    tabLegendre

Table of values of normalized Legendre associated functions
See `computePlmcostheta` of the `SphericalHarmonics` to see the normalization.
"""
const tabLegendre = zeros(Float64,N+1,nbK)


# tables which depend on x,q

"""
    tabOmega_serial

Table of values of the angular velocity of the system.
"""
tabOmega_serial = zeros(Float64,nbK)

"""
    tabKappaSqOverTwoOmega_serial

Table of values of the ratio of `kappa^2/(2 Omega)`.
"""
tabKappaSqOverTwoOmega_serial = zeros(Float64,nbK)

"""
    tabAln_serial

Table of values of the matrix `Aln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabAln_serial = zeros(Float64,N+1,N+1)

"""
    tabBln_serial

Table of values of the matrix `Bln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabBln_serial = zeros(Float64,N+1,N+1)

"""
    tabCln_serial

Table of values of the matrix `Cln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabCln_serial = zeros(Float64,N+1,N+1)

"""
    tabDln_serial

Table of values of the matrix `Dln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabDln_serial = zeros(Float64,N+1,N+1)

"""
    tabFln_serial

Table of values of the matrix `Fln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabFln_serial = zeros(Float64,N+1,N+1)

"""
    tabGln_serial

Table of values of the matrix `Gln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabGln_serial = zeros(Float64,N+1,N+1)

"""
    tabHln_serial

Table of values of the matrix `Hln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
tabHln_serial = zeros(Float64,N+1,N+1)

###############################################
# Fill Iln table
###############################################

"""
    I_seed!()

Initializes `tabIln[1,1]` with the value `I(m,m,m)`.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
function I_seed!()
    prod = 1/(3/4+m+1)
    for k=1:m
        prod *= 2*(2*k+1)/(3/4+m+1+k)
    end
    tabIln[1,1] = prod
end

"""
    I_firstLine!()

Initializes the first column `tabIln[k,1]` with the values `I(m,l,m)` for k=2,...,2N+1.
Initializes the first lign `tabIln[1,k]` by symmetry.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
function I_firstLine!()
    current = tabIln[1,1]
    for l=m+1:m+2*N
        current *= (l-3/4-m-1)/(l+3/4+m+1)*sqrt(((l+m)*(2*l+1))/((l-m)*(2*l-1)))
        tabIln[l-m+1,1] = current
        tabIln[1,l-m+1] = current #symmetry
   end
end

"""
    I_fill!()

Computes the table `tabIln` by recursion.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
function I_fill!()
    for n=m+1:m+N
        for l=n:m+2*N-(n-m)
            Ip = (sqrt(((l+m+1)*(l-m+1))/((2*l+1)*(2*l+3)))*tabIln[l-m+1+1,n-m+1-1]
                 + sqrt(((l+m)*(l-m))/((2*l+1)*(2*l-1))) *tabIln[l-m+1-1,n-m+1-1])
            tabIln[l-m+1,n-m+1] = sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)))*Ip
            if (n>=m+2)
                tabIln[l-m+1,n-m+1] -= sqrt(((n+m-1)*(n-m-1)*(2*n+1))/((n+m)*(n-m)*(2*n-3)))*tabIln[l-m+1,n-m+1-2]
            end
            tabIln[n-m+1,l-m+1] = tabIln[l-m+1,n-m+1] #symmetry
        end
    end
end

"""
    tabIln!()

Computes the whole table `tabIln`.
See `Aoki et al. (1979)` (a=3/4).
"""
function tabIln!()
    I_seed!()
    I_firstLine!()
    I_fill!()
end

"""
    tabIln_clear!()

Clears the whole table `tabIln`.
"""
function tabIln_clear!()
    for i=1:2*N+1
        for j=1:2*N+1
            tabIln[i,j] = 0.0
        end
    end
end

###############################################
# Fill Ialm table : Ia(m-1;l,m-1) with l=m-1,...,m+2N-1
###############################################

"""
    Ia_seed!()

Initializes `tabIln[1,1]` with the value `Ia(m-1,m-1,m-1)`.
See `Aoki et al. (1979)` for a definition (a=7/4).
"""
function Ia_seed!()
    prod = 1/(7/4+m-1+1)
    for k=1:m-1
        prod *= 2*(2*k+1)/(7/4+m-1+1+k)
    end
    tabIaln[1] = prod
end

"""
    Ia_firstLine!()

Initializes the table `tabIaln[k]` with the values `Ia(m-1,l,m-1)` for k=2,...,2N+1.
See `Aoki et al. (1979)` for a definition (a=7/4).
"""
function Ia_firstLine!()
    current = tabIaln[1]
    for l=m:m+2*N-1
        current *= (l-7/4-m+1-1)/(l+7/4+m-1+1)*sqrt(((l+m-1)*(2*l+1))/((l-m+1)*(2*l-1)))
        tabIaln[l-m+2] = current
   end
end

"""
    tabIaln!()

Computes the whole table `tabIaln`.
See `Aoki et al. (1979)` (a=7/4).
"""
function tabIaln!()
    Ia_seed!()
    Ia_firstLine!()
end

"""
    tabIaln_clear!()

Clears the whole table `tabIaln`.
"""
function tabIaln_clear!()
    for i=1:2*N+1
        tabIaln[i] = 0.0
    end
end

###############################################
# Fill Jln table
###############################################

"""
    J_seed!()

Initializes `tabJln[1,1]` with the value `J(m,m,m)`.
Initializes `tabJln[2,1]` with the value `J(m,m+1,m)` and `tabJln[1,2]` by symmetry.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
function J_seed!()
     tabJln[1,1] = tabIln[1,1]*(3/4+2*m+1)/(m)
     tabJln[2,1] = -(7/4)*sqrt(2*m+3)/(3/4+2*m+1)*tabJln[1,1]
     tabJln[1,2] = tabJln[2,1] #symmetry
end

"""
    I_firstLine!()

Initializes the first column `tabJln[k,1]` with the values `J(m,l,m)` for k=3,...,2N+1.
Initializes the first lign `tabJln[1,k]` by symmetry.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
function J_firstLine!()
    for l=m+2:m+2*N
        current = (sqrt(((l-m)*(l-m-1)*(2*l+1))/((l+m)*(l+m-1)*(2*l-3)))*tabJln[l-m+1-2,1]
                  + 4*sqrt(((2*m+1)*(2*l+1)*(2*l-1))/((2*m)*(l+m)*(l+m-1)))*tabIaln[l-m+2-1])
        tabJln[l-m+1,1] = current
        tabJln[1,l-m+1] = current #symmetry
   end
end

"""
    J_fill!()

Computes the table `tabJln` by recursion.
See `Aoki et al. (1979)` for a definition (a=3/4).
"""
function J_fill!()
    for n=m+1:m+N
        for l=n:m+2*N-(n-m)
            Jp = (sqrt(((l+m+1)*(l-m+1))/((2*l+1)*(2*l+3)))*tabJln[l-m+1+1,n-m+1-1]
                 + sqrt(((l+m)*(l-m))/((2*l+1)*(2*l-1))) *tabJln[l-m+1-1,n-m+1-1])
            tabJln[l-m+1,n-m+1] = sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)))*Jp
            if (n>=m+2)
                tabJln[l-m+1,n-m+1] -= sqrt(((n+m-1)*(n-m-1)*(2*n+1))/((n+m)*(n-m)*(2*n-3)))*tabJln[l-m+1,n-m+1-2]
            end
            tabJln[n-m+1,l-m+1] = tabJln[l-m+1,n-m+1] #symmetry
        end
    end
end

"""
    tabJln!()

Computes the whole table `tabJln`.
See `Aoki et al. (1979)` (a=3/4).
"""
function tabJln!()
    J_seed!()
    J_firstLine!()
    J_fill!()
end

"""
    tabJln_clear!()

Clears the whole table `tabJln`.
"""
function tabJln_clear!()
    for i=1:2*N+1
        for j=1:2*N+1
            tabJln[i,j] = 0.0
        end
    end
end

###############################################
# Fill tabLegendre table
###############################################

"""
    tabLegendre!()

Fills the table `tabLegendre` with mid-point rule evaluation of the normalized Legendre associated functions.
See `computePlmcostheta` of the `SphericalHarmonics` for a definition.
"""
function tabLegendre!()
    for k=1:nbK
        xik = -1 +(2/nbK)*(k-0.5)
        thk = acos(xik)
        P = computePlmcostheta(thk,m+N)
        for l=m:m+N
            tabLegendre[l-m+1,k] = sqrt(pi)*P[(l,m)]
        end
    end
end

"""
    tabLegendre_clear!()

Clears the table `tabLegendre`.
"""
function tabLegendre_clear!()
    for i=1:N+1
        for k=1:nbK
            tabLegendre[i,k] = 0.0
        end
    end
end

###############################################
# Fill tabOmega table
###############################################

"""
    tabOmega!(x,q,[args])

Fills the table `tabOmega` with mid-point rule evaluation of the dimensionless angular velocity.
"""
function tabOmega!(x::Float64=0.0, q::Float64=1.0,tabOmega=tabOmega_serial)
    for k=1:nbK
        xik = -1 +(2/nbK)*(k-0.5)
        tabOmega[k] = Omega(xik,x,q)
    end
end

"""
    tabOmega_clear!([args])

Clears the table `tabOmega`.
"""
function tabOmega_clear!(tabOmega=tabOmega_serial)
    for k=1:nbK
        tabOmega[k] = 0.0
    end
end

###############################################
# Fill tabAlphaSqOverTwoOmega table
###############################################

"""
    tabKappaSqOverTwoOmega!()

Fills the table `tabOmega` with mid-point rule evaluation of the dimensionless ratio `kappa^2/(2Omega)`.
"""
function tabKappaSqOverTwoOmega!(x::Float64=0.0, q::Float64=1.0,tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial)
    for k=1:nbK
        xik = -1 +(2/nbK)*(k-0.5)
        tabKappaSqOverTwoOmega[k] = kappaSqOverTwoOmega(xik,x,q)
    end
end

"""
    tabKappaSqOverTwoOmega_clear!([args])

Clears the table `tabKappaSqOverTwoOmega`.
"""
function tabKappaSqOverTwoOmega_clear!(tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial)
    for k=1:nbK
        tabKappaSqOverTwoOmega[k] = 0.0
    end
end


###############################################
# Fill matrix elements table
###############################################

"""
    tabAlnFln!([args])

Fills the tables `tabAln` and `tabFln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
function tabAlnFln!(tabAln=tabAln_serial,tabFln=tabFln_serial,tabOmega=tabOmega_serial)
    for l=m:m+N
        for n=l:m+N
            integral = 0
            for k=1:nbK
                integral += tabLegendre[l-m+1,k]*tabOmega[k]*tabLegendre[n-m+1,k]
            end
            integral *= 2/nbK
            tabAln[l-m+1,n-m+1] = abs(m)*integral
            tabFln[l-m+1,n-m+1] = 2*integral
            if (n>l)
                tabAln[n-m+1,l-m+1] = tabAln[l-m+1,n-m+1]
                tabFln[n-m+1,l-m+1] = tabFln[l-m+1,n-m+1]
            end
        end
    end
end

"""
    tabBln!(x,[args])

Fills the tables `tabBln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
function tabBln!(x::Float64=0.0,tabBln=tabBln_serial)
    for l=m:m+N
        for n=m:m+N
            tabBln[l-m+1,n-m+1] = sqrt(1-x)*((1/2)*(sqrt(((2*l+1)*(l+m+1)*(l-m+1))/(2*l+3))*tabJln[l-m+1+1,n-m+1]
                                + tabJln[l-m+1,n-m+1] ))
            if (l>m)
                tabBln[l-m+1,n-m+1] -= sqrt(1-x)*(1/2)*sqrt(((2*l+1)*(l+m)*(l-m))/(2*l-1))*tabJln[l-m+1-1,n-m+1]
            end
        end
    end
end

"""
    tabCln!(x,[args])

Fills the tables `tabCln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
function tabCln!(x::Float64=0.0,tabCln=tabCln_serial)
    for l=m:m+N
        for n=m:m+N
            tabCln[l-m+1,n-m+1] = m*sqrt(1-x)*tabJln[l-m+1,n-m+1]
        end
    end
end

"""
    tabDln!(x,[args])

Fills the tables `tabDln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
function tabDln!(x::Float64=0.0, tabDln=tabDln_serial)
    for l=m:m+N
        for n=m:m+N
            tabDln[l-m+1,n-m+1] = sqrt(1-x)*((1/2)*(1/(2*n+1)-(eps/3)/((1-x)^(2/3)))*((-sqrt(((2*n+1)*(n+m+1)*(n-m+1))/(2*n+3))
                                *tabIln[l-m+1,n-m+1+1])- tabIln[l-m+1,n-m+1] ))
            if (n>m)
                tabDln[l-m+1,n-m+1] += sqrt(1-x)*((1/2)*(1/(2*n+1)-(eps/3)/((1-x)^(2/3)))* (sqrt(((2*n+1)*(n+m)*(n-m))/(2*n-1))
                                       *tabIln[l-m+1,n-m+1-1]))
            end
        end
    end
end

"""
    tabGln!(x,[args])

Fills the tables `tabGln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
function tabGln!(x::Float64=0.0,tabGln=tabGln_serial)
    for l=m:m+N
        for n=m:m+N
            tabGln[l-m+1,n-m+1] = -m*sqrt(1-x)*(1/(2*n+1)-(eps/3)/((1-x)^(2/3)))*tabIln[l-m+1,n-m+1]
        end
    end
end

"""
    tabHln!([args])

Fills the tables `tabHln`.
See `Aoki et al. (1979)` for the case `x=0` and `q=1`.
"""
function tabHln!(tabHln=tabHln_serial,tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial)
    for l=m:m+N
        for n=l:m+N
            integral = 0
            for k=1:nbK
                integral += tabLegendre[l-m+1,k]*tabKappaSqOverTwoOmega[k]*tabLegendre[n-m+1,k]
            end
            integral *= 2/nbK
            tabHln[l-m+1,n-m+1] = integral
            if (n>l)
                tabHln[n-m+1,l-m+1] = tabHln[l-m+1,n-m+1]
            end
        end
    end
end


###############################################
# Fill matrix elements table
###############################################

"""
    table_constant_fill!()

Fills the recursion tables of `I`, `Ia` and `J`.
"""
function table_constant_fill!()
    tabIln!()
    tabIaln!()
    tabJln!()
    tabLegendre!()
end

"""
    table_function_fill!(x,q,[args])

Fills the function tables (angular velocity and epicyclic frequency).
"""
function table_function_fill!(x::Float64=0.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial)
    tabOmega!(x,q,tabOmega)
    tabKappaSqOverTwoOmega!(x,q,tabKappaSqOverTwoOmega)
end

"""
    matrix_fill!(x,q,[args])

Fills the matrix elements.
"""
function matrix_fill!(x::Float64=0.0, q::Float64=1.0,
    tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
    tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
    tabHln=tabHln_serial, tabOmega=tabOmega_serial,
    tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial)
    tabAlnFln!(tabAln,tabFln,tabOmega)
    tabBln!(x,tabBln)
    tabCln!(x,tabCln)
    tabDln!(x,tabDln)
    tabGln!(x,tabGln)
    tabHln!(tabHln,tabKappaSqOverTwoOmega)
end

######################################
# Clear table
######################################

"""
    tabAln_clear!([args])

Clears the matrix elements of `Aln`.
"""
function tabAln_clear!(tabAln=tabAln_serial)
    for i=1:N+1
        for j=1:N+1
            tabAln[i,j] = 0.0
        end
    end
end

"""
    tabBln_clear!([args])

Clears the matrix elements of `Bln`.
"""
function tabBln_clear!(tabBln=tabBln_serial)
    for i=1:N+1
        for j=1:N+1
            tabBln[i,j] = 0.0
        end
    end
end

"""
    tabCln_clear!([args])

Clears the matrix elements of `Cln`.
"""
function tabCln_clear!(tabCln=tabCln_serial)
    for i=1:N+1
        for j=1:N+1
            tabCln[i,j] = 0.0
        end
    end
end

"""
    tabDln_clear!([args])

Clears the matrix elements of `Dln`.
"""
function tabDln_clear!(tabDln=tabDln_serial)
    for i=1:N+1
        for j=1:N+1
            tabDln[i,j] = 0.0
        end
    end
end

"""
    tabFln_clear!([args])

Clears the matrix elements of `Fln`.
"""
function tabFln_clear!(tabFln=tabFln_serial)
    for i=1:N+1
        for j=1:N+1
            tabFln[i,j] = 0.0
        end
    end
end

"""
    tabGln_clear!([args])

Clears the matrix elements of `Gln`.
"""
function tabGln_clear!(tabGln=tabGln_serial)
    for i=1:N+1
        for j=1:N+1
            tabGln[i,j] = 0.0
        end
    end
end

"""
    tabHln_clear!([args])

Clears the matrix elements of `Hln`.
"""
function tabHln_clear!(tabHln=tabHln_serial)
    for i=1:N+1
        for j=1:N+1
            tabHln[i,j] = 0.0
        end
    end
end

"""
    table_function_clear!([args])

Clears the function tables.
"""
function table_function_clear!(tabOmega=tabOmega_serial,
    tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial)
    tabOmega_clear!(tabOmega)
    tabKappaSqOverTwoOmega_clear!(tabKappaSqOverTwoOmega)
end

"""
    matrix_clear!([args])

Clears all the matrix tables.
"""
function matrix_clear!(tabAln=tabAln_serial, tabBln=tabBln_serial,
    tabCln=tabCln_serial, tabDln=tabDln_serial, tabFln=tabFln_serial,
    tabGln=tabGln_serial, tabHln=tabHln_serial)
    tabAln_clear!(tabAln)
    tabBln_clear!(tabBln)
    tabCln_clear!(tabCln)
    tabDln_clear!(tabDln)
    tabFln_clear!(tabFln)
    tabGln_clear!(tabGln)
    tabHln_clear!(tabHln)
end

######################################
# Initialize constant Tables
######################################

table_constant_fill!()
