///   *** Matrix Exponentiation [2^3*log(n)]
///   by calling exponent this will calculate (ara^p)%md;
///   here n is the dimention of the matrix ara

#define ll long long
#define DIM 5

ll tmp[DIM][DIM],result[DIM][DIM];

void cpy(ll ar1[][DIM],ll ar2[][DIM],ll n){
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            ar1[i][j] = ar2[i][j];
}

void mult(ll ar1[][DIM], ll ar2[][DIM], ll n, ll mod){
    int i,j,k;
    memset(tmp,0,sizeof tmp);
    for(i=1;i<=n;i++){
        for(j=1;j<=n;j++){
            for(k=1;k<=n;k++){
                tmp[i][j] = (tmp[i][j] + ar1[i][k]*ar2[k][j])%mod;
            }
        }
    }
    cpy(ar1,tmp,n);
}

void exponent(ll ara[][DIM], ll p,ll n, ll mod){
    int i,j;
    //initializing identity
    for(i=1;i<=n;i++)
        for(j=1;j<=n;j++)
            if(i==j) result[i][j] = 1;
            else result[i][j] = 0;

    while(p){
        if(p&1){
            mult(result,ara,n,mod);
        }
        mult(ara,ara,n,mod);
        p >>= 1;
    }
    cpy(ara,result,n);
}

