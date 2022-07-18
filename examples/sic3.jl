#

using QP
using Oscar

N = 3; R = QQX(N); Rp = QQXp(N); I = Im(N) + Ih(N); I2 = Im(N) + Ih2(N); I0 = Itr0(N); trx = TrX(N); trx2 = TrX2(N);  Rt = QQXt(N); #S = ProjectiveScheme(R,I);

Rs, ph = quo(R,Im(N)); 