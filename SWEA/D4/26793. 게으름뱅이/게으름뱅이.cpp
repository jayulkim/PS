#include <bits/stdc++.h>
#define all(v) v.begin(), v.end()
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
typedef long double ld;
typedef tuple<ll, ll, ll> lll;
typedef tuple<ll, ll, ll, ll> llll;
typedef tuple<ll, ll, ll, ll, ll> lllll;
typedef tuple<ll, ll, ll, ll, ll, ll> llllll;
typedef pair<ll, ll> pll;
typedef pair<string, ll> psl;
ll n = 0, m = 0, k = 0;
class Fastmod {
private:
	using ll = long long;
	using ull = unsigned long long;
	ull mod;
	ull norm(ll x) {
		ll m = (ll)mod;
		x %= m;
		if (x < 0) x += m;
		return (ull)x;
	}
	ull addmod_u(ull a, ull b) {
		if (a >= mod - b) {
			return a - (mod - b);
		}
		return a + b;
	}
public:
	Fastmod(ull mod) : mod(mod) {}
	ull addmod(ull a, ull b) {
		return addmod_u(a, b);
	}
	ull addmod(ll a, ll b) {
		return addmod_u(norm(a), norm(b));
	}
	ull addmod(const vector<ull>& v) {
		ull result = 0;
		for (auto& i : v) {
			result = addmod_u(result, i);
		}
		return result;
	}
	ll addmod(const vector<ll>& v) {
		ll result = 0;
		for (auto& i : v) {
			result = addmod(result, i);
		}
		return result;
	}
	ull fastmulmod(ull a, ull b) {
		ll result = a * b - mod * (ull)(1.L / mod * a * b);
		if (result < 0) {
			result += mod;
		}
		if ((ull)result >= mod) {
			result -= mod;
		}
		return (ull)result;
	}
	ull fastmulmod(ll a, ll b) {
		return fastmulmod(norm(a), norm(b));
	}
	ull fastmulmod(const vector<ull>& v) {
		ull result = 1;
		for (auto& i : v) {
			result = fastmulmod(result, i);
		}
		return result;
	}
	ull fastmulmod(const vector<ll>& v) {
		ull result = 1;
		for (auto& i : v) {
			result = fastmulmod((ll)result, i);
		}
		return result;
	}
	ull mulmod(ull a, ull b) {
		ull result = 0;
		while (b) {
			if (b & 1) {
				result = addmod_u(result, a);
			}
			a = addmod_u(a, a);
			b >>= 1;
		}
		return result;
	}
	ull mulmod(ll a, ll b) {
		return mulmod(norm(a), norm(b));
	}
	ull mulmod(const vector<ull>& v) {
		ull result = 1;
		for (auto& i : v) {
			result = mulmod(result, i);
		}
		return result;
	}
	ull mulmod(const vector<ll>& v) {
		ull result = 1;
		for (auto& i : v) {
			result = mulmod((ll)result, i);
		}
		return result;
	}
	ull fastpowmod(ll a, ull b) {
		ull x = norm(a);
		ull result = 1;
		while (b) {
			if (b & 1) {
				result = fastmulmod(result, x);
			}
			x = fastmulmod(x, x);
			b >>= 1;
		}
		return result;
	}
	ull powmod(ll a, ull b) {
		ull x = norm(a);
		ull result = 1;
		while (b) {
			if (b & 1) {
				result = mulmod(result, x);
			}
			x = mulmod(x, x);
			b >>= 1;
		}
		return result;
	}
	ull divmod(ll up, ll down) {
		return fastmulmod(norm(up), fastpowmod(norm(down), mod - 2));
	}
};
class Vectormath {
private:
	using ll = long long;
	using ull = unsigned long long;
public:
	ull add(const vector<ull>& v) {
		ull result = 0;
		for (auto& i : v) {
			result += i;
		}
		return result;
	}
	ll add(const vector<ll>& v) {
		ll result = 0;
		for (auto& i : v) {
			result += i;
		}
		return result;
	}
	ull mul(const vector<ull>& v) {
		ull result = 1;
		for (auto& i : v) {
			result *= i;
		}
		return result;
	}
	ll mul(const vector<ll>& v) {
		ll result = 1;
		for (auto& i : v) {
			result *= i;
		}
		return result;
	}
};
class Matrix {
private:
	using ll = long long;
	using ull = unsigned long long;
	using matrix = vector<vector<ll>>;
	using matrix_u = vector<vector<ull>>;
	ull mod;
	Fastmod fm;
public:
	Matrix(ull mod) : mod(mod), fm(mod) {}
	ll pow(ll a, ull b) {
		ll result = 1;
		while (b) {
			if (b & 1) {
				result *= a;
			}
			a *= a;
			b >>= 1;
		}
		return result;
	}
	ull pow(ull a, ull b) {
		ull result = 1;
		while (b) {
			if (b & 1) {
				result *= a;
			}
			a *= a;
			b >>= 1;
		}
		return result;
	}
	matrix addmod(const matrix& a, const matrix& b) {
		ll n = a.size();
		ll m = a[0].size();
		if (n != (ll)b.size() || m != (ll)b[0].size()) {
			return { {} };
		}
		matrix result = a;
		for (ll i = 0; i < n; i++) {
			for (ll j = 0; j < m; j++) {
				result[i][j] = fm.addmod(result[i][j], b[i][j]);
			}
		}
		return result;
	}
	matrix_u addmod(const matrix_u& a, const matrix_u& b) {
		ll n = a.size();
		ll m = a[0].size();
		if (n != (ll)b.size() || m != (ll)b[0].size()) {
			return { {} };
		}
		matrix_u result = a;
		for (ll i = 0; i < n; i++) {
			for (ll j = 0; j < m; j++) {
				result[i][j] = fm.addmod(result[i][j], b[i][j]);
			}
		}
		return result;
	}
	matrix mulmod(const matrix& a, const matrix& b) {
		if ((ll)a[0].size() != (ll)b.size()) {
			return { {} };
		}
		matrix result(a.size(), vector<ll>(b[0].size(), 0));
		for (ll i = 0; i < (ll)a.size(); i++) {
			for (ll k = 0; k < (ll)a[0].size(); k++) {
				if (a[i][k]) {
					for (ll j = 0; j < (ll)b[0].size(); j++) {
						result[i][j] = fm.addmod(result[i][j], (ll)fm.fastmulmod(a[i][k], b[k][j]));
					}
				}
			}
		}
		return result;
	}
	matrix mul(const matrix& a, const matrix& b) {
		if ((ll)a[0].size() != (ll)b.size()) {
			return { {} };
		}
		matrix result(a.size(), vector<ll>(b[0].size(), 0));
		for (ll i = 0; i < (ll)a.size(); i++) {
			for (ll k = 0; k < (ll)a[0].size(); k++) {
				if (a[i][k]) {
					for (ll j = 0; j < (ll)b[0].size(); j++) {
						result[i][j] += a[i][k] * b[k][j];
					}
				}
			}
		}
		return result;
	}
	matrix_u mulmod(const matrix_u& a, const matrix_u& b) {
		if (a[0].size() != b.size()) {
			return { {} };
		}
		matrix_u result(a.size(), vector<ull>(b[0].size(), 0));
		for (ll i = 0; i < (ll)a.size(); i++) {
			for (ll k = 0; k < (ll)a[0].size(); k++) {
				if (a[i][k]) {
					for (ll j = 0; j < (ll)b[0].size(); j++) {
						result[i][j] = fm.addmod(result[i][j], fm.fastmulmod(a[i][k], b[k][j]));
					}
				}
			}
		}
		return result;
	}
	matrix_u mul(const matrix_u& a, const matrix_u& b) {
		if (a[0].size() != b.size()) {
			return { {} };
		}
		matrix_u result(a.size(), vector<ull>(b[0].size(), 0));
		for (ll i = 0; i < (ll)a.size(); i++) {
			for (ll k = 0; k < (ll)a[0].size(); k++) {
				if (a[i][k]) {
					for (ll j = 0; j < (ll)b[0].size(); j++) {
						result[i][j] += a[i][k] * b[k][j];
					}
				}
			}
		}
		return result;
	}
	matrix identity(ll n) {
		matrix result(n, vector<ll>(n, 0));
		for (ll i = 0; i < n; i++) {
			result[i][i] = 1;
		}
		return result;
	}
	matrix powmod(matrix a, ull b) {
		if (a.size() != a[0].size()) {
			return { {} };
		}
		matrix result = identity((ll)a.size());
		while (b) {
			if (b & 1) {
				result = mulmod(result, a);
			}
			a = mulmod(a, a);
			b >>= 1;
		}
		return result;
	}
	matrix pow(matrix a, ull b) {
		if (a.size() != a[0].size()) {
			return { {} };
		}
		matrix result = identity((ll)a.size());
		while (b) {
			if (b & 1) {
				result = mul(result, a);
			}
			a = mul(a, a);
			b >>= 1;
		}
		return result;
	}
	matrix_u identity(ull n) {
		matrix_u result(n, vector<ull>(n, 0));
		for (ull i = 0; i < n; i++) {
			result[i][i] = 1;
		}
		return result;
	}
	matrix_u powmod(matrix_u a, ull b) {
		if (a.size() != a[0].size()) {
			return { {} };
		}
		matrix_u result = identity((ull)a.size());
		while (b) {
			if (b & 1) {
				result = mulmod(result, a);
			}
			a = mulmod(a, a);
			b >>= 1;
		}
		return result;
	}
	matrix_u pow(matrix_u a, ull b) {
		if (a.size() != a[0].size()) {
			return { {} };
		}
		matrix_u result = identity((ull)a.size());
		while (b) {
			if (b & 1) {
				result = mul(result, a);
			}
			a = mul(a, a);
			b >>= 1;
		}
		return result;
	}
};
//class Primecheck {
//private:
//	using ull = unsigned long long;
//public:
//	bool bigisprime(ull x) {
//		static const ull sample[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
//		if (x < 2 || x % 6 % 4 != 1) {
//			return (x | 1) == 3;
//		}
//		Fastmod fm(x);
//		ull s = __builtin_ctzll(x - 1);
//		ull d = (x - 1) >> s;
//		for (ull a : sample) {
//			ull temp = fm.fastpowmod(a % x, d);
//			ull i = s;
//			while (temp != 1 && temp != x - 1 && a % x && i--) {
//				temp = fm.fastmulmod(temp, temp);
//			}
//			if (temp != x - 1 && i != s) {
//				return false;
//			}
//		}
//		return true;
//	}
//	bool isprime(ull x) {
//		if (x < 2) {
//			return false;
//		}
//		if (x % 2 == 0) {
//			return x == 2;
//		}
//		for (ull i = 3; i <= x / i; i += 2) {
//			if (x % i == 0) {
//				return false;
//			}
//		}
//		return true;
//	}
//	vector<bool> eratosthenes(ull n) {
//		vector<bool>result(n + 1, false);
//		result[1] = true;
//		for (int i = 2; i * i <= n; i++) {
//			if (!result[i]) {
//				for (int j = i * i; j <= n; j += i) {
//					result[j] = true;
//				}
//			}
//		}
//		return result;
//	}
//};
//class Fastgcd {
//private:
//	using ll = long long;
//	using ull = unsigned long long;
//	using lll = tuple<ll, ll, ll>;
//public:
//	ull fastgcd(ull a, ull b) {
//		ull shift = 0;
//		if (a == 0) {
//			return b;
//		}
//		if (b == 0) {
//			return a;
//		}
//		shift = __builtin_ctzll(a | b);
//		a >>= __builtin_ctzll(a);
//		do {
//			b >>= __builtin_ctzll(b);
//			if (a > b) {
//				swap(a, b);
//			}
//			b -= a;
//		} while (b);
//		return a << shift;
//	}
//	ull fastgcd(const vector<ull>& v) {
//		ull result = 0;
//		for (auto& i : v) {
//			result = fastgcd(result, i);
//		}
//		return result;
//	}
//	ull fastgcd(const vector<ll>& v) {
//		ull result = 0;
//		for (auto& i : v) {
//			result = fastgcd(result, i);
//		}
//		return result;
//	}
//	ull fastlcm(ull a, ull b) {
//		if (!a) {
//			return b;
//		}
//		if (!b) {
//			return a;
//		}
//		return (a * b) / fastgcd(a, b);
//	}
//	ull fastlcm(const vector<ull>& v) {
//		ull result = 0;
//		for (auto& i : v) {
//			result = fastlcm(result, i);
//		}
//		return result;
//	}
//	ull fastlcm(const vector<ll>& v) {
//		ull result = 0;
//		for (auto& i : v) {
//			result = fastlcm(result, i);
//		}
//		return result;
//	}
//	lll eea(ll a, ll mod) {
//		ll r1 = a, r2 = mod;
//		ll s1 = 1, s2 = 0;
//		ll t1 = 0, t2 = 1;
//		while (r2) {
//			ll q = r1 / r2;
//			ll r = r1 - q * r2;
//			r1 = r2;
//			r2 = r;
//			ll s = s1 - q * s2;
//			s1 = s2;
//			s2 = s;
//			ll t = t1 - q * t2;
//			t1 = t2;
//			t2 = t;
//		}
//		return { r1, s1, t1 };
//	}
//};
//class Primefactor {
//private:
//	using ll = long long;
//	using ull = unsigned long long;
//	using pull = pair<ull, ull>;
//	using pll = pair<ll, ll>;
//	Fastgcd fg;
//	Primecheck pc;
//public:
//	ull pollardrho(ull x) {
//		if (x % 2 == 0) {
//			return 2;
//		}
//		Fastmod fm(x);
//		auto f = [&](ull a) {
//			return (fm.fastmulmod(a, a) + 1) % x;
//			};
//		ull a = 0, b = 0;
//		ull t = 30;
//		ull prd = 2;
//		ull i = 1;
//		ull q = 0;
//		while (t++ % 40 || fg.fastgcd(prd, x) == 1) {
//			if (a == b) {
//				a = ++i;
//				b = f(a);
//			}
//			q = fm.fastmulmod(prd, (a > b ? a - b : b - a));
//			if (q) {
//				prd = q;
//			}
//			a = f(a);
//			b = f(f(b));
//			if (!prd) {
//				return x;
//			}
//		}
//		return fg.fastgcd(prd, x);
//	}
//	vector<ull> bigfactor(ull x) {
//		vector<ull> result;
//		vector<ull> v = { x };
//		while (!v.empty()) {
//			ull temp = v.back();
//			v.pop_back();
//			if (temp == 1) {
//				continue;
//			}
//			if (pc.bigisprime(temp)) {
//				result.push_back(temp);
//				continue;
//			}
//			ull d = pollardrho(temp);
//			v.push_back(d);
//			v.push_back(temp / d);
//		}
//		return result;
//	}
//	void bigfactor(ull x, vector<ull>& result) {
//		vector<ull> v = { x };
//		while (!v.empty()) {
//			ull temp = v.back();
//			v.pop_back();
//			if (temp == 1) {
//				continue;
//			}
//			if (pc.bigisprime(temp)) {
//				result.push_back(temp);
//				continue;
//			}
//			ull d = pollardrho(temp);
//			v.push_back(d);
//			v.push_back(temp / d);
//		}
//	}
//	vector<ll> bigfactor(ll x) {
//		vector<ll> result;
//		vector<ll> v = { x };
//		while (!v.empty()) {
//			ll temp = v.back();
//			v.pop_back();
//			if (temp == 1) {
//				continue;
//			}
//			if (pc.bigisprime(temp)) {
//				result.push_back(temp);
//				continue;
//			}
//			ll d = pollardrho(temp);
//			v.push_back(d);
//			v.push_back(temp / d);
//		}
//		return result;
//	}
//	void bigfactor(ll x, vector<ll>& result) {
//		vector<ll> v = { x };
//		while (!v.empty()) {
//			ll temp = v.back();
//			v.pop_back();
//			if (temp == 1) {
//				continue;
//			}
//			if (pc.bigisprime(temp)) {
//				result.push_back(temp);
//				continue;
//			}
//			ll d = pollardrho(temp);
//			v.push_back(d);
//			v.push_back(temp / d);
//		}
//	}
//	vector<ull> factor(ull x) {
//		vector<ull> result;
//		for (ull i = 2; i <= x / i; i++) {
//			while (x % i == 0) {
//				result.push_back(i);
//				x /= i;
//			}
//		}
//		if (x > 1) {
//			result.push_back(x);
//		}
//		return result;
//	}
//	void factor(ull x, vector<ull>& result) {
//		for (ull i = 2; i <= x / i; i++) {
//			while (x % i == 0) {
//				result.push_back(i);
//				x /= i;
//			}
//		}
//		if (x > 1) {
//			result.push_back(x);
//		}
//	}
//	vector<ll> factor(ll x) {
//		vector<ll> result;
//		for (ull i = 2; i <= x / i; i++) {
//			while (x % i == 0) {
//				result.push_back(i);
//				x /= i;
//			}
//		}
//		if (x > 1) {
//			result.push_back(x);
//		}
//		return result;
//	}
//	void factor(ll x, vector<ll>& result) {
//		for (ll i = 2; i <= x / i; i++) {
//			while (x % i == 0) {
//				result.push_back(i);
//				x /= i;
//			}
//		}
//		if (x > 1) {
//			result.push_back(x);
//		}
//	}
//	vector<pll> factorialfactor(ll n, const vector<ll>& prime) {
//		vector<pll> result;
//		for (auto& i : prime) {
//			ll count = 0;
//			for (ll j = i; j <= n;) {
//				count += n / j;
//				if (j > n / i) {
//					break;
//				}
//				j *= i;
//			}
//			if (count > 0) {
//				result.push_back({ i, count });
//			}
//		}
//		return result;
//	}
//	void factorialfactor(ll n, const vector<ll>& prime, vector<ll>& result) {
//		for (auto& i : prime) {
//			ll count = 0;
//			for (ll j = i; j <= n;) {
//				count += n / j;
//				if (j > n / i) {
//					break;
//				}
//				j *= i;
//			}
//			if (count > 0) {
//				result[i] += count;
//			}
//		}
//	}
//	vector<ull> getbigdivisor(ull n) {
//		vector<ull> f = bigfactor(n);
//		sort(f.begin(), f.end());
//		vector<pull> pf;
//		for (auto& i : f) {
//			if (pf.empty() || pf.back().first != i) {
//				pf.push_back({ i, 1 });
//			}
//			else {
//				pf.back().second++;
//			}
//		}
//		vector<ull> result = { 1 };
//		for (auto& [a, b] : pf) {
//			size_t size = result.size();
//			ull temp = 1;
//			for (int i = 1; i <= b; i++) {
//				temp *= a;
//				for (size_t j = 0; j < size; j++) {
//					result.push_back(result[j] * temp);
//				}
//			}
//		}
//		return result;
//	}
//	void getbigdivisor(ull n, vector<ull>& result) {
//		vector<ull> f = bigfactor(n);
//		sort(f.begin(), f.end());
//		vector<pull> pf;
//		for (auto& i : f) {
//			if (pf.empty() || pf.back().first != i) {
//				pf.push_back({ i, 1 });
//			}
//			else {
//				pf.back().second++;
//			}
//		}
//		vector<ull> v = { 1 };
//		for (auto& [a, b] : pf) {
//			size_t size = v.size();
//			ull temp = 1;
//			for (int i = 1; i <= b; i++) {
//				temp *= a;
//				for (size_t j = 0; j < size; j++) {
//					v.push_back(v[j] * temp);
//				}
//			}
//		}
//		for (auto& i : v) {
//			result.push_back(i);
//		}
//	}
//	vector<ll> getbigdivisor(ll n) {
//		vector<ll> f = bigfactor(n);
//		sort(f.begin(), f.end());
//		vector<pll> pf;
//		for (auto& i : f) {
//			if (pf.empty() || pf.back().first != i) {
//				pf.push_back({ i, 1 });
//			}
//			else {
//				pf.back().second++;
//			}
//		}
//		vector<ll> result = { 1 };
//		for (auto& [a, b] : pf) {
//			size_t size = result.size();
//			ll temp = 1;
//			for (int i = 1; i <= b; i++) {
//				temp *= a;
//				for (size_t j = 0; j < size; j++) {
//					result.push_back(result[j] * temp);
//				}
//			}
//		}
//		return result;
//	}
//	void getbigdivisor(ll n, vector<ll>& result) {
//		vector<ll> f = bigfactor(n);
//		sort(f.begin(), f.end());
//		vector<pll> pf;
//		for (auto& i : f) {
//			if (pf.empty() || pf.back().first != i) {
//				pf.push_back({ i, 1 });
//			}
//			else {
//				pf.back().second++;
//			}
//		}
//		vector<ll> v = { 1 };
//		for (auto& [a, b] : pf) {
//			size_t size = v.size();
//			ll temp = 1;
//			for (int i = 1; i <= b; i++) {
//				temp *= a;
//				for (size_t j = 0; j < size; j++) {
//					v.push_back(v[j] * temp);
//				}
//			}
//		}
//		for (auto& i : v) {
//			result.push_back(i);
//		}
//	}
//	vector<ull> getdivisor(ull n) {
//		vector<ull> result;
//		for (ull i = 1; i <= n / i; i++) {
//			if (n % i == 0) {
//				result.push_back(i);
//				if (i != n / i) {
//					result.push_back(n / i);
//				}
//			}
//		}
//		return result;
//	}
//	void getdivisor(ull n, vector<ull>& result) {
//		for (ull i = 1; i <= n / i; i++) {
//			if (n % i == 0) {
//				result.push_back(i);
//				if (i != n / i) {
//					result.push_back(n / i);
//				}
//			}
//		}
//	}
//	vector<ll> getdivisor(ll n) {
//		vector<ll> result;
//		for (ll i = 1; i <= n / i; i++) {
//			if (n % i == 0) {
//				result.push_back(i);
//				if (i != n / i) {
//					result.push_back(n / i);
//				}
//			}
//		}
//		return result;
//	}
//	void getdivisor(ll n, vector<ll>& result) {
//		for (ll i = 1; i <= n / i; i++) {
//			if (n % i == 0) {
//				result.push_back(i);
//				if (i != n / i) {
//					result.push_back(n / i);
//				}
//			}
//		}
//	}
//};
//class Eulerphi {
//private:
//	using ull = unsigned long long;
//	Primefactor pf;
//public:
//	ull eulerphi(ull n) {
//		ull result = 1;
//		for (ull i = 2; i <= n / i; i++) {
//			if (n % i == 0) {
//				while (n % i == 0) {
//					result *= i;
//					n /= i;
//				}
//				result /= i;
//				result *= (i - 1);
//			}
//		}
//		if (n > 1) {
//			result *= (n - 1);
//		}
//		return result;
//	}
//	ull bigeulerphi(ull n) {
//		vector<ull> fac = pf.bigfactor(n);
//		sort(fac.begin(), fac.end());
//		fac.erase(unique(fac.begin(), fac.end()), fac.end());
//		ull result = n;
//		for (ull p : fac) {
//			result = result / p * (p - 1);
//		}
//		return result;
//	}
//};
class Combination {
private:
	using ull = unsigned long long;
	ull mod;
	Fastmod fm;
public:
	Combination(ull mod) : mod(mod), fm(mod) {}
	ull nCr(ull n, ull r) {
		ull u = 1, d = 1;
		if (r > n) {
			return 0;
		}
		if (r == 0 || n == r) {
			return 1;
		}
		r = min(r, n - r);
		for (int i = 0; i < r; i++) {
			u = fm.fastmulmod(u, n);
			n--;
		}
		for (ull i = 1; i <= r; i++) {
			d = fm.fastmulmod(d, i);
		}
		return fm.fastmulmod(u, fm.fastpowmod(d, mod - 2));
	}
	ull nHr(ull n, ull r) {
		if (n == 0 && r > 0) {
			return 0;
		}
		return nCr(n + r - 1, r);
	}
	ull fastnCr(ull n, ull r, const vector<ull>& factorial, const vector<ull>& invfactorial) {
		if (r > n) {
			return 0;
		}
		if (r == 0 || n == r) {
			return 1;
		}
		return fm.fastmulmod(factorial[n], fm.fastmulmod(invfactorial[r], invfactorial[n - r]));
	}
	ull fastnHr(ull n, ull r, const vector<ull>& factorial, const vector<ull>& invfactorial) {
		if (n == 0 && r > 0) {
			return 0;
		}
		return fastnCr(n + r - 1, r, factorial, invfactorial);
	}
	ull lucasnCr(ull n, ull r) {
		if (r > n) {
			return 0;
		}
		ull result = 1;
		while (n > 0 || r > 0) {
			ull ni = n % mod;
			ull ri = r % mod;
			if (ri > ni) {
				return 0;
			}
			result = fm.fastmulmod(result, nCr(ni, ri));
			n /= mod;
			r /= mod;
		}
		return result;
	}
	ull lucasnHr(ull n, ull r) {
		if (n == 0 && r > 0) {
			return 0;
		}
		return lucasnCr(n + r - 1, r);
	}
};
class Dijkstra {
private:
	using ll = long long;
	using pll = pair<ll, ll>;
	ll n;
	vector<vector<pll>>graph;
public:
	Dijkstra(ll n, const vector<vector<pll>>& graph) : n(n), graph(graph) {}
	vector<ll> run(const vector<ll>& start) {
		priority_queue<pll, vector<pll>, greater<pll>>pq;
		vector<ll>dist(n + 1, LLONG_MAX);
		for (auto& i : start) {
			dist[i] = 0;
			pq.push({ 0, i });
		}
		while (!pq.empty()) {
			auto [a, b] = pq.top();
			pq.pop();
			if (dist[b] >= a) {
				for (auto& i : graph[b]) {
					if (dist[i.first] > dist[b] + i.second) {
						dist[i.first] = dist[b] + i.second;
						pq.push({ dist[i.first], i.first });
					}
				}
			}
		}
		return dist;
	}
	pair<vector<ll>, vector<vector<pll>>> revrun(const vector<ll>& start) {
		priority_queue<pll, vector<pll>, greater<pll>>pq;
		vector<vector<pll>>revgraph(n + 1);
		vector<ll>dist(n + 1, LLONG_MAX);
		for (auto& i : start) {
			dist[i] = 0;
			pq.push({ 0, i });
		}
		while (!pq.empty()) {
			auto [a, b] = pq.top();
			pq.pop();
			if (dist[b] >= a) {
				for (auto& i : graph[b]) {
					if (dist[i.first] > dist[b] + i.second) {
						vector<pll>().swap(revgraph[i.first]);
						revgraph[i.first].push_back({ b, i.second });
						dist[i.first] = dist[b] + i.second;
						pq.push({ dist[i.first], i.first });
					}
					else if (dist[i.first] == dist[b] + i.second) {
						revgraph[i.first].push_back({ b, i.second });
					}
				}
			}
		}
		return { dist, revgraph };
	}
};
class Toposort {
private:
	using ll = long long;
	ll n;
	vector<vector<ll>>dag;
public:
	Toposort(ll n, const vector<vector<ll>>& dag) : n(n), dag(dag) {}
	pair<bool, vector<ll>> toposort() {
		vector<ll>degree(n + 1, 0);
		for (int i = 1; i <= n; i++) {
			for (auto& j : dag[i]) {
				degree[j]++;
			}
		}
		vector<ll>result;
		queue<ll>q;
		for (int i = 1; i <= n; i++) {
			if (!degree[i]) {
				q.push(i);
			}
		}
		while (!q.empty()) {
			ll temp = q.front();
			q.pop();
			result.push_back(temp);
			for (auto& i : dag[temp]) {
				if (--degree[i] == 0) {
					q.push(i);
				}
			}
		}
		return { result.size() == n, result };
	}
};
class Dsu {
private:
	using ll = long long;
	using pll = pair<ll, ll>;
	ll n;
	vector<ll>parent;
	vector<ll>rank;
public:
	void init() {
		parent.resize(n + 1);
		iota(parent.begin(), parent.end(), 0);
		rank.assign(n + 1, 1);
	}
	Dsu(ll n) : n(n) {
		init();
	}
	ll getparent(ll x) {
		if (x == parent[x]) {
			return x;
		}
		return parent[x] = getparent(parent[x]);
	}
	void merge(ll a, ll b) {
		a = getparent(a);
		b = getparent(b);
		if (a == b) {
			return;
		}
		if (rank[a] > rank[b]) {
			rank[a] += rank[b];
			rank[b] = 0;
			parent[b] = a;
		}
		else {
			rank[b] += rank[a];
			rank[a] = 0;
			parent[a] = b;
		}
	}
	bool hascycle(const vector<vector<ll>>& graph) {
		init();
		for (int i = 1; i <= n; i++) {
			for (auto& j : graph[i]) {
				if (i == j) {
					return true;
				}
				if (i < j) {
					if (getparent(i) == getparent(j)) {
						return true;
					}
					merge(i, j);
				}
			}
		}
		return false;
	}
	bool hascycle(const vector<vector<pll>>& graph) {
		init();
		for (int i = 1; i <= n; i++) {
			for (auto& j : graph[i]) {
				if (i == j.first) {
					return true;
				}
				if (i < j.first) {
					if (getparent(i) == getparent(j.first)) {
						return true;
					}
					merge(i, j.first);
				}
			}
		}
		return false;
	}
};
class Mst {
private:
	using ll = long long;
	using lll = tuple<ll, ll, ll>;
	ll n;
public:
	Mst(ll n) : n(n) {}
	ll minkruskal(priority_queue<lll, vector<lll>, greater<lll>> pq) {
		Dsu ds(n);
		ll result = 0;
		ll count = 0;
		while (!pq.empty() && count < n - 1) {
			auto [a, b, c] = pq.top();
			pq.pop();
			if (ds.getparent(b) != ds.getparent(c)) {
				count++;
				result += a;
				ds.merge(b, c);
			}
		}
		return (count == n - 1 ? result : LLONG_MAX);
	}
	ll maxkruskal(priority_queue<lll> pq) {
		Dsu ds(n);
		ll result = 0;
		ll count = 0;
		while (!pq.empty() && count < n - 1) {
			auto [a, b, c] = pq.top();
			pq.pop();
			if (ds.getparent(b) != ds.getparent(c)) {
				count++;
				result += a;
				ds.merge(b, c);
			}
		}
		return (count == n - 1 ? result : LLONG_MAX);
	}
};
class Floydwarshall {
private:
	using ll = long long;
	ll n;
	vector<vector<ll>>graph;
public:
	Floydwarshall(ll n, const vector<vector<ll>>& graph) : n(n), graph(graph) {}
	vector<vector<ll>> run() {
		vector<vector<ll>>dist = graph;
		for (int i = 1; i <= n; i++) {
			dist[i][i] = 0;
		}
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				for (int k = 1; k <= n; k++) {
					if (dist[j][i] != LLONG_MAX && dist[i][k] != LLONG_MAX) {
						dist[j][k] = min(dist[j][k], dist[j][i] + dist[i][k]);
					}
				}
			}
		}
		return dist;
	}
	pair<vector<vector<ll>>, vector<vector<ll>>>revrun() {
		vector<vector<ll>>dist = graph;
		vector<vector<ll>>next(n + 1, vector<ll>(n + 1, -1));
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				if (graph[i][j] != LLONG_MAX) {
					next[i][j] = j;
				}
			}
		}
		for (int i = 1; i <= n; i++) {
			dist[i][i] = 0;
			next[i][i] = i;
		}
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				for (int k = 1; k <= n; k++) {
					if (dist[j][i] != LLONG_MAX && dist[i][k] != LLONG_MAX) {
						ll temp = dist[j][i] + dist[i][k];
						if (temp < dist[j][k]) {
							dist[j][k] = temp;
							next[j][k] = next[j][i];
						}
					}
				}
			}
		}
		return { dist, next };
	}
	vector<ll> getpath(ll start, ll end, const vector<vector<ll>>& next) {
		if (next[start][end] == -1) {
			return {};
		}
		vector<ll>path = { start };
		for (int i = 0; i <= n; i++) {
			if (start == end) {
				return path;
			}
			start = next[start][end];
			path.push_back(start);
		}
		return {};
	}
};
int main(void) {
	ios::sync_with_stdio(false);
	cin.tie(0);
	ll num = 0;
	cin >> num;
	while (num--) {
		cin >> n;
		map<ll, ll>Map;
		for (int i = 0; i < n; i++) {
			ll a = 0, b = 0;
			cin >> a >> b;
			Map[b] += a;
		}
		ll Min = LLONG_MAX;
		ll temp = 0;
		for (auto& [a, b] : Map) {
			temp += b;
			Min = min(Min, a - temp);
		}
		cout << Min << '\n';
	}
	return 0;
}