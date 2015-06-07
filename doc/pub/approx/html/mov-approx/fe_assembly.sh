#!/bin/sh
# Create one HTML file with different types of assembly illustrations
out=fe_assembly.html
cat > $out <<EOF
<html>
<head>
</head>
<body>
EOF
dirs="fe_assembly_regular_2x2 fe_assembly_regular_4x4 fe_assembly_irregular"
for dir in $dirs; do
name=`echo $dir | sed -e 's/fe_assembly_//g'`
cat >> $out <<EOF

<h1>Mesh with $name numbering</h1>

<img src="$dir/fe_mesh.png" width=600>
<p>
EOF

cp $dir/index.html .
doconce replace "\"frame_" "\"$dir/frame_" index.html
doconce grab --from '<script language="Javascript">' --to- '</body>' index.html >> $out
done

cat >> $out <<EOF
</body>
</html>
EOF
