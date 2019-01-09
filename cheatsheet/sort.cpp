std::sort(duplets.begin(), duplets.end(),
          [](duplet<T> x, duplet<T> y) { return x.ind < y.ind; });