// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Template specialization of model parameters for the model described in "twoband_Cu.hpp".

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::FeAs5Orb<PointGroup>>> {
public:
  ModelParameters() = default;

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  auto get_U() const {
    return U_;
  }
  auto get_V() const {
    return V_;
  }
  auto get_J() const {
    return J_;
  }
  auto get_Jp() const {
    return Jp_;
  }

private:

  double U_ = 0;
  double V_ = 0;

  double J_ = 0;
  double Jp_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::FeAs5Orb<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(J_);
  buffer_size += concurrency.get_buffer_size(Jp_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::FeAs5Orb<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, J_);
  concurrency.pack(buffer, buffer_size, position, Jp_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::FeAs5Orb<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, J_);
  concurrency.unpack(buffer, buffer_size, position, Jp_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::FeAs5Orb<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("FeAs5Orb-model");

    auto try_to_read = [&](const std::string& name, auto& val) {
      try {
        reader_or_writer.execute(name, val);
      }
      catch (const std::exception& r_e) {
      }
    };

    try_to_read("U", U_);
    try_to_read("V", V_);

    try_to_read("J", J_);
    try_to_read("Jp", Jp_);

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    throw(std::logic_error("No model parameters in input file."));
  }
}
