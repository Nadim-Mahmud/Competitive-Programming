#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<ll> vll;
#define ntm(i, n) for(int i = 0; i<(n); ++i)
#define rep(i, a, b) for(int i = a; i<(b); ++i)
#define all(x) (x).begin(),(x).end()
#define trav(x, v) for(auto &x : (v))
#define sz(x) (int)(x).size()
#define mp(a, b) make_pair(a, b)
#define X first
#define Y second
#define pb(x) push_back(x)
#define PS(x) cout << " " << (x)
#define P1(x) cout << (x) << endl
#define P2(x,y) cout << (x) << " " << (y) << endl
#define P3(x,y,z) cout << (x) << " " << (y) << " " << (z) << endl
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define PF(x, y) fixed << setprecision(y) << x

using namespace std;

struct TwoSat {
    int N;
    vector<vi> gr;
    vi values; // 0 = false , 1 = true
    TwoSat(int n = 0) : N(n), gr(2*n) {}
    int add_var() {
        gr.emplace_back();
        gr.emplace_back();
        return N++;
    }
    void either(int f, int j) {
        f = max(2*f, -1-2*f);
        j = max(2*j, -1-2*j);
        //P3("either: ", f, j);
        gr[f^1].push_back(j);
        gr[j^1].push_back(f);
    }
    void set_value(int x) { either(x, x); }
    void at_most_one(const vi& li) {
        if (sz(li) <= 1) return;
        int cur = ~li[0];
        rep(i,2,sz(li)) {
            int next = add_var();
            either(cur, ~li[i]);
            either(cur, next);
            either(~li[i], next);
            cur = ~next;
        }
        either(cur, ~li[1]);
    }
    vi val, comp, z; int time = 0;
    int dfs(int i) {
        int low = val[i] = ++time, x; z.push_back(i);
        trav(e, gr[i]) if (!comp[e])
            low = min(low, val[e] ?: dfs(e));
        ++time;
        if (low == val[i]) do {
            x = z.back(); z.pop_back();
            comp[x] = time;
            if (values[x>>1] == -1)
                values[x>>1] = !(x&1);
        } while (x != i);
        return val[i] = low;
    }
    bool solve() {
        values.assign(N, -1);
        val.assign(2*N, 0); comp = val;
        rep(i,0,2*N) if (!comp[i]) dfs(i);
        rep(i,0,N) if (comp[2*i] == comp[2*i+1]) return 0;
        return 1;
    }
};
