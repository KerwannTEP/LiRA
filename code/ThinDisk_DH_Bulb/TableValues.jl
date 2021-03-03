using SphericalHarmonics
using StaticArrays

# Tables shifted by 1

# tables which don't depend on x,q
const tabIln = zeros(Float64,2*N+1,2*N+1)
const tabIaln = zeros(Float64,2*N+1)
const tabJln = zeros(Float64,2*N+1,2*N+1)
const tabLegendre = zeros(Float64,N+1,nbK)

# tables which depend on x,q
tabOmega = zeros(Float64,nbK)
tabAlphaSqOverTwoOmega = zeros(Float64,nbK)

tabAln = zeros(Float64,N+1,N+1)
tabBln = zeros(Float64,N+1,N+1)
tabCln = zeros(Float64,N+1,N+1)
tabDln = zeros(Float64,N+1,N+1)
tabFln = zeros(Float64,N+1,N+1)
tabGln = zeros(Float64,N+1,N+1)
tabHln = zeros(Float64,N+1,N+1)

# apply in correct order

###############################################
# Fill Iln table
###############################################

# initialize I(m,m)
function I_seed!()
    prod = 1/(3/4+m+1)
    for k=1:m
        prod *= 2*(2*k+1)/(3/4+m+1+k)
    end
    tabIln[1,1] = prod
end

# l=m+1,...,m+2N
function I_firstLine!()
    current = tabIln[1,1]
    for l=m+1:m+2*N
        current *= (l-3/4-m-1)/(l+3/4+m+1)*sqrt(((l+m)*(2*l+1))/((l-m)*(2*l-1)))
        tabIln[l-m+1,1] = current
        tabIln[1,l-m+1] = current #symmetry
   end
end

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

function tabIln!()
    I_seed!()
    I_firstLine!()
    I_fill!()
end

function tabIln_clear!()
    for i=1:2*N+1
        for j=1:2*N+1
            tabIln[i,j] = 0.0
        end
    end
end

###############################################
# Fill Ialm table : Ia(l,m-1) with l=m-1,...,m+2N-1
###############################################

# initialize Ia(m-1,m-1) (a=7/4)
function Ia_seed!()
    prod = 1/(7/4+m-1+1)
    for k=1:m
        prod *= 2*(2*k+1)/(7/4+m-1+1+k)
    end
    tabIaln[1] = prod
end

# l=m,...,m+2N-1
function Ia_firstLine!()
    current = tabIaln[1]
    for l=m:m+2*N-1
        current *= (l-7/4-m-1-1)/(l+7/4+m-1+1)*sqrt(((l+m-1)*(2*l+1))/((l-m+1)*(2*l-1)))
        tabIaln[l-m+2] = current
   end
end

function tabIaln!()
    Ia_seed!()
    Ia_firstLine!()
end

function tabIaln_clear!()
    for i=1:2*N+1
        tabIaln[i] = 0.0
    end
end

###############################################
# Fill Jln table
###############################################

# initialize J(m,m) and J(m+1,m)
function J_seed!()
    prod = 1
    for k=1:m
        prod *= 2*(2*k+1)/(3/4+m+k)
    end
    tabJln[1,1] = prod/m
    prod *= -(7/4)*sqrt(2*m+3)/(3/4+2*m+1)
    tabJln[2,1] = prod
    tabJln[2,1] = tabJln[2,1] #symmetry
end

# l=m+1,...,m+2N
function J_firstLine!()
    for l=m+2:m+2*N
        current = (sqrt(((l-m)*(l-m-1)*(2*l+1))/((l+m)*(l+m-1)*(2*l-3)))*tabJln[l-m+1-2,1]
                  + sqrt(((2*m+1)*(2*l+1)*(2*l-1))/((2*m)*(l+m)*(l+m-1)))*tabIaln[l-m+2-1])
        tabJln[l-m+1,1] = current
        tabJln[1,l-m+1] = current #symmetry
   end
end

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

function tabJln!()
    J_seed!()
    J_firstLine!()
    J_fill!()
end

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

function tabLegendre!()
    for k=1:nbK
        xik = -1 +(2/nbK)*(k-0.5)
        thk = acos(xik)
        P = computePlmcostheta(thk,m+N,m)
        for l=m:m+N
            tabLegendre[l-m+1,k] = sqrt(pi)*P[(l,m)]
        end
    end
end

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

function tabOmega!(x::Float64=1.0, q::Float64=1.0)
    for k=1:nbK
        xik = -1 +(2/nbK)*(k-0.5)
        tabOmega[k] = Omega(xik,x,q)
    end
end

function tabOmega_clear!()
    for k=1:nbK
        tabOmega[k] = 0.0
    end
end

###############################################
# Fill tabAlphaSqOverTwoOmega table
###############################################

function tabAlphaSqOverTwoOmega!(x::Float64=1.0, q::Float64=1.0)
    for k=1:nbK
        xik = -1 +(2/nbK)*(k-0.5)
        tabAlphaSqOverTwoOmega[k] = alphaSqOverTwoOmega(xik,x,q)
    end
end

function tabAlphaSqOverTwoOmega_clear!()
    for k=1:nbK
        tabAlphaSqOverTwoOmega[k] = 0.0
    end
end


###############################################
# Fill matrix elements table
###############################################

function tabAlnFln!()
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

function tabBln!()
    for l=m:m+N
        for n=l:m+N
            tabBln[l-m+1,n-m+1] = (1/2)*(sqrt(((2*l+1)*(l+m+1)*(l-m+1))/(2*l+3))*tabJln[l-m+1+1,n-m+1]
                                + tabJln[l-m+1,n-m+1] )
            if (l>m)
                tabBln[l-m+1,n-m+1] -= (1/2)*sqrt(((2*l+1)*(l+m)*(l-m))/(2*l-1))*tabJln[l-m+1-1,n-m+1]
            end
            if (n>l)
                tabBln[n-m+1,l-m+1] = tabBln[l-m+1,n-m+1]
            end
        end
    end
end

function tabCln!()
    for l=m:m+N
        for n=l:m+N
            tabCln[l-m+1,n-m+1] = m*tabJln[l-m+1,n-m+1]
            if (n>l)
                tabCln[n-m+1,l-m+1] = tabCln[l-m+1,n-m+1]
            end
        end
    end
end


function tabDln!(x::Float64=1.0, q::Float64=1.0)
    for l=m:m+N
        for n=l:m+N
            tabDln[l-m+1,n-m+1] = (1/2)*(1/(2*n+1)-(eps/3)*(q/x)^(2/3))*(-sqrt(((2*n+1)*(n+m+1)*(n-m+1))/(2*n+3))*tabIln[l-m+1+1,n-m+1]
                                - tabIln[l-m+1,n-m+1] )
            if (l>m)
                tabDln[l-m+1,n-m+1] += (1/2)*(1/(2*n+1)-(eps/3)*(q/x)^(2/3))* sqrt(((2*n+1)*(n+m)*(n-m))/(2*n-1))*tabIln[l-m+1-1,n-m+1]
            end
            if (n>l)
                tabDln[n-m+1,l-m+1] = tabDln[l-m+1,n-m+1]
            end
        end
    end
end

function tabGln!(x::Float64=1.0, q::Float64=1.0)
    for l=m:m+N
        for n=l:m+N
            tabGln[l-m+1,n-m+1] = (1/2)*(1/(2*n+1)-(eps/3)*(q/x)^(2/3))*tabIln[l-m+1,n-m+1]
            if (n>l)
                tabGln[n-m+1,l-m+1] = tabGln[l-m+1,n-m+1]
            end
        end
    end
end

function tabHln!()
    for l=m:m+N
        for n=l:m+N
            integral = 0
            for k=1:nbK
                integral += tabLegendre[l-m+1,k]*tabAlphaSqOverTwoOmega[k]*tabLegendre[n-m+1,k]
            end
            integral *= 2/nbK
            tabHln[l-m+1,n-m+1] = integral
            if (n>l)
                tabHln[n-m+1,l-m+1] = tabHln[l-m+1,n-m+1]
            end
        end
    end
end

# do in order

function table_constant_fill!()
    tabIln!()
    tabIaln!()
    tabJln!()
    tabLegendre!()
end

function table_function_fill!(x::Float64=1.0, q::Float64=1.0)
    tabOmega!(x,q)
    tabAlphaSqOverTwoOmega!(x,q)
end

function matrix_fill!(x::Float64=1.0, q::Float64=1.0)
    tabAlnFln!()
    tabBln!()
    tabCln!()
    tabDln!(x,q)
    tabGln!(x,q)
    tabHln!()
end

######################################
# Clear table
######################################

function tabAln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabAln[i,j] = 0.0
        end
    end
end

function tabBln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabBln[i,j] = 0.0
        end
    end
end

function tabCln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabCln[i,j] = 0.0
        end
    end
end

function tabDln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabDln[i,j] = 0.0
        end
    end
end

function tabFln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabFln[i,j] = 0.0
        end
    end
end

function tabGln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabGln[i,j] = 0.0
        end
    end
end

function tabHln_clear!()
    for i=1:N+1
        for j=1:N+1
            tabHln[i,j] = 0.0
        end
    end
end

function table_function_clear!()
    tabOmega_clear!()
    tabAlphaSqOverTwoOmega_clear!()
end

function matrix_clear!()
    tabAln_clear!()
    tabBln_clear!()
    tabCln_clear!()
    tabDln_clear!()
    tabFln_clear!()
    tabGln_clear!()
    tabHln_clear!()
end

######################################
# Initialize constant Tables
######################################

table_constant_fill!()
