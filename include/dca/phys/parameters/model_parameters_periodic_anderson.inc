// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Model parameters for the two-dimensional periodic Anderson model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::PeriodicAndersonModel<PointGroup>>> {
public:
  ModelParameters() : t_(0.), epsilon_f_(0.), V_(0.), V_prime_(0.), U_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t() const {
    return t_;
  }
  void set_t(const double t) {
    t_ = t;
  }

  double get_epsilon_f() const {
    return epsilon_f_;
  }
  void set_epsilon_f(const double epsilon_f) {
    epsilon_f_ = epsilon_f;
  }

  double get_V() const {
    return V_;
  }
  void set_V(const double V) {
    V_ = V;
  }

  double get_V_prime() const {
    return V_prime_;
  }
  void set_V_prime(const double V_prime) {
    V_prime_ = V_prime;
  }

  double get_U() const {
    return U_;
  }
  void set_U(const double U) {
    U_ = U;
  }

private:
  double t_;
  double epsilon_f_;
  double V_;
  double V_prime_;
  double U_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::PeriodicAndersonModel<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  return concurrency.get_buffer_size(t_) + concurrency.get_buffer_size(epsilon_f_) +
         concurrency.get_buffer_size(V_) + concurrency.get_buffer_size(V_prime_) +
         concurrency.get_buffer_size(U_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::PeriodicAndersonModel<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, const int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, epsilon_f_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
  concurrency.pack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::PeriodicAndersonModel<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, const int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, epsilon_f_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
  concurrency.unpack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::PeriodicAndersonModel<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  reader_or_writer.open_group("periodic-Anderson-model");

  try {
    reader_or_writer.execute("t", t_);
  }
  catch (const std::exception&) {
  }
  try {
    reader_or_writer.execute("epsilon-f", epsilon_f_);
  }
  catch (const std::exception&) {
  }
  try {
    reader_or_writer.execute("V", V_);
  }
  catch (const std::exception&) {
  }
  try {
    reader_or_writer.execute("V-prime", V_prime_);
  }
  catch (const std::exception&) {
  }
  try {
    reader_or_writer.execute("U", U_);
  }
  catch (const std::exception&) {
  }

  reader_or_writer.close_group();

  if (ReaderOrWriter::is_reader) {
    std::cout << "\n\n"
              << "\tperiodic-Anderson-model :\n"
              << "\t-----------------------\n\n"
              << "\t\t t         : " << t_ << "\n"
              << "\t\t epsilon-f : " << epsilon_f_ << "\n"
              << "\t\t V         : " << V_ << "\n"
              << "\t\t V-prime   : " << V_prime_ << "\n"
              << "\t\t U         : " << U_ << "\n\n"
              << std::endl;
  }
}
