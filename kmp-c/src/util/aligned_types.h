#pragma once
#include <vector>
#include "aligned_allocator.h"

template<typename T, Alignment U = AlignAVX> using aligned_vector = std::vector<T, aligned_allocator<T, U>>;
