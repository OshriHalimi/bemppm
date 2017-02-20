#ifndef FMM_COMMON_HPP
#define FMM_COMMON_HPP

#include "../common/types.hpp"
#include "../common/shared_ptr.hpp"

namespace fmm {

template <typename T> using Vector = Bempp::Vector<T>;
template <typename T> using Matrix = Bempp::Matrix<T>;
template <typename T> using Point3D = Bempp::Point3D<T>;

using Bempp::shared_ptr;
using Bempp::make_shared_from_ref;
using Bempp::make_shared_from_const_ref;
using Bempp::null_deleter;
using Bempp::dynamic_pointer_cast;
using Bempp::static_pointer_cast;
using Bempp::const_pointer_cast;

}

#endif
