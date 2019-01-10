static sparse_vec conv(const sparse_vec &a, const sparse_vec &b)
{
    sparse_vec out(a.len + b.len - 1);
    for (auto x : a.duplets)
    {
        for (auto y : b.duplets)
        {
            out.append(x.ind + y.ind, x.val * y.val);
        }
    }
    return out;
}