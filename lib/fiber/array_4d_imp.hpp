#include <stdexcept>

namespace Fiber
{

template <typename T>
inline Array4D<T>::Array4D()
{
    m_extents[0] = 0;
    m_extents[1] = 0;
    m_extents[2] = 0;
    m_extents[3] = 0;
    m_storage = 0;
    m_owns = false;
}

template <typename T>
inline Array4D<T>::Array4D(int extent0, int extent1, int extent2, int extent3)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2, extent3);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_extents[3] = extent3;
    m_storage = new T[extent0 * extent1 * extent2 * extent3];
    m_owns = true;
}

template <typename T>
inline Array4D<T>::Array4D(int extent0, int extent1, int extent2, int extent3, T* data)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2, extent3);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_extents[3] = extent3;
    m_storage = data;
    m_owns = false;
}

template <typename T>
inline Array4D<T>::~Array4D()
{
    if (m_owns && m_storage)
        delete[] m_storage;
}

template <typename T>
inline T& Array4D<T>::operator()(int index0, int index1, int index2, int index3)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1, index2, index3);
#endif
    return m_storage[
            index0 +
            m_extents[0] * index1 +
            m_extents[0] * m_extents[1] * index2 +
            m_extents[0] * m_extents[1] * m_extents[2] * index3];
}

template <typename T>
inline const T& Array4D<T>::operator()(int index0, int index1, int index2, int index3) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1, index2, index3);
#endif
    return m_storage[
            index0 +
            m_extents[0] * index1 +
            m_extents[0] * m_extents[1] * index2 +
            m_extents[0] * m_extents[1] * m_extents[2] * index3];
}

template <typename T>
inline int Array4D<T>::extent(int dimension) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_dimension(dimension);
#endif
    return m_extents[dimension];
}

template <typename T>
inline void Array4D<T>::set_size(int extent0, int extent1, int extent2, int extent3)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2, extent3);
#endif
    if (extent0 * extent1 * extent2 * extent3 ==
            m_extents[0] * m_extents[1] * m_extents[2] * m_extents[3]) {
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
        m_extents[3] = extent3;
    }
    else {
        if (m_owns && m_storage) {
            delete[] m_storage;
            m_storage = 0;
        }
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
        m_extents[3] = extent3;
        m_storage = new T[extent0 * extent1 * extent2 * extent3];
        m_owns = true;
    }
}

template <typename T>
inline typename Array4D<T>::iterator Array4D<T>::begin()
{
    return m_storage;
}

template <typename T>
inline typename Array4D<T>::const_iterator Array4D<T>::begin() const
{
    return m_storage;
}

template <typename T>
inline typename Array4D<T>::iterator Array4D<T>::end()
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2] * m_extents[3];
}

template <typename T>
inline typename Array4D<T>::const_iterator Array4D<T>::end() const
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2] * m_extents[3];
}

#ifdef FIBER_CHECK_ARRAY_BOUNDS
template <typename T>
inline void Array4D<T>::check_dimension(int dimension) const
{
    if (dimension < 0 || 3 < dimension)
        throw std::invalid_argument("Invalid dimension");
}

template <typename T>
inline void Array4D<T>::check_extents(int extent0, int extent1, int extent2, int extent3) const
{
    if (extent0 <= 0 || extent1 <= 0 || extent2 <= 0 || extent3 <= 0)
        throw std::length_error("Invalid extent");
}

template <typename T>
inline void Array4D<T>::check_indices(int index0, int index1, int index2, int index3) const
{
    if (index0 < 0 || m_extents[0] <= index0 ||
        index1 < 0 || m_extents[1] <= index1 ||
        index2 < 0 || m_extents[2] <= index2 ||
        index3 < 0 || m_extents[3] <= index3)
        throw std::out_of_range("Invalid index");
}
#endif // FIBER_CHECK_ARRAY_BOUNDS

} // namespace Fiber
