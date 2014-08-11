#pragma once
#include <vector>
#include "aligned_allocator.h"

template<typename T> using aligned_vector = std::vector<T, aligned_allocator<T>>;
