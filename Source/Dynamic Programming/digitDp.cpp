/**
*	Description : generate all number upto the number give into the string st[]
				Given implementation will count how many numbers and those numbers digits sum divisible by n
*   Complexity : O(multiplicatins of parameters)
*	Use : set n = mod number | set ln = length(st)
*/


char st[MX];
int n,ln;
int dp[11][3][101][101];

int bct(int pos,bool lst,int sum,int mult){
    if(pos>=ln){
        if(sum%n==0&&mult%n==0) return 1;
        else return 0;
    }
    if(dp[pos][lst][sum][mult]!=-1)
        return dp[pos][lst][sum][mult];
    int x=0,y;
    for(int i=0; i<=9; i++){
        if(lst){
            y = st[pos]-'0';
            if(i>y) break;
            x += bct(pos+1,i==y,(sum+i)%n,(mult*10+i)%n);
        }
        else{
            x += bct(pos+1,0,(sum+i)%n,(mult*10+i)%n);
        }
    }
    return dp[pos][lst][sum][mult] = x;
}
