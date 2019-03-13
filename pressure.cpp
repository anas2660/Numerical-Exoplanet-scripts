#include "exoterDE.cpp"
int main() {
  ExoEq::LoadData();
  vector<thread*> threads = {
      new thread(ExoterDE, "CoRoT7b", 4.8, 0.8, 1.68, 0.09),
      new thread(ExoterDE, "GJ1214b", 6.55, 0.98, 2.678, 0.13)};
  for (auto thr : threads) thr->join();
  if (PLOT) system("python ternaryplot.py CoRoT7b GJ1214b");
}
