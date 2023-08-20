function vn = Replicate(v,n)

vn = repmat( v, n, 1 );
vn = vn(:)';

return;