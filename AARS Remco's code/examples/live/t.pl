$s = "R_01 R_02 R_03 R_04 R_05 R_06 R_07 R_08 R_09 R_10 R_11 R_12 R_13 R_14 C_01 C_02 C_03 C_05 C_06 C_08 C_09 C_10 C_11 C_12 C_13 C_14 E_01 E_02 E_03 E_04 E_05 E_06 E_07 E_08 E_09 E_10 E_11 E_12 E_13 E_14 I_01 I_02 I_03 I_04 I_05 I_06 I_07 I_08 I_09 I_10 I_11 I_12 I_13 I_14 L_01 L_02 L_03 L_04 L_05 L_06 L_07 L_08 L_09 L_10 L_11 L_12 L_13 L_14 K_01 K_02 K_03 K_04 K_05 K_07 K_08 K_10 K_11 K_13 K_14 M_01 M_02 M_03 M_04 M_05 M_06 M_07 M_08 M_09 M_10 M_11 M_12 M_13 M_14 W_01 W_02 W_22 W_03 W_04 W_05 W_07 W_09 W_10 W_11 W_12 W_13 W_14 W_Na Y_01 Y_02 Y_03 Y_04 Y_05 Y_06 Y_07 Y_08 Y_09 Y_19 Y_10 Y_11 Y_12 Y_13 Y_14 V_01 V_02 V_03 V_04 V_05 V_06 V_07 V_08 V_09 V_10 V_11 V_12 V_13 V_14 A_01 A_02 A_03 A_04 A_05 A_06 A_07 A_09 A_10 A_12 A_13 A_14 N_09 N_10 N_11 N_13 N_14 D_01 D_02 D_03 D_04 D_05 D_06 D_07 D_08 D_09 D_10 D_11 D_12 D_13 D_14 G_01 G_02 G_03 G_04 G_05 G_06 G_07 G_08 G_09 G_10 G_11 G_12 G_13 G_14 H_01 H_02 H_03 H_04 H_05 H_06 H_07 H_08 H_09 H_10 H_11 H_12 H_13 H_14 F_01 F_02 F_03 F_05 F_06 F_08 F_10 F_12 F_13 F_14 P_01 P_02 P_03 P_04 P_05 P_06 P_07 P_08 P_09 P_10 P_11 P_12 P_13 P_14 S_01 S_02 S_05 S_06 S_09 S_10 S_11 S_12 S_13 S_14 T_01 T_02 T_03 T_04 T_05 T_06 T_07 T_08 T_09 T_10 T_11 T_12 T_13 T_14";
@s =split(" ", $s);



$x = 'VGADPESLRTIQNCYHMWKF';
for ($i = 0; $i < 20; $i++) {
    $order{substr($x,$i,1)} = $i;
}
$pre = 'R';
for $s (@s) {
    $p = substr($s,0,1);
    if($pre ne $p) {
        print "\n        ";
        $pre = $p;
    }
    print "$s=$order{$p},"
}
exit;


print '             <distribution id="R.prior" monophyletic="true" spec="beast.math.distributions.MRCAPrior" tree="@Tree.t:AARS_20_1">
                <taxonset id="R" spec="TaxonSet">
';
$pre = 'R';
for $s (@s) {
    $p = substr($s,0,1);
    if($pre ne $p) {
print'                </taxonset>
            </distribution>
';
print '             <distribution id="'.$p.'.prior" monophyletic="true" spec="beast.math.distributions.MRCAPrior" tree="@Tree.t:AARS_20_1">
                <taxonset id="'.$p.'" spec="TaxonSet">
';
        $pre = $p;
    }
    print "                   <taxon id=\"$s\"/>\n";
}
print'                </taxonset>
            </distribution>
';





exit;
while ($s = <>) {
    if ($s =~ /\t([^\t]+)\t([^\t]+)/) {
        $s = "\t".substr($1,0,6)."\t$2";
        $s =~ s/Ala/A/;
        $s =~ s/Arg/R/;
        $s =~ s/Asn/N/;
        $s =~ s/Asp/D/;
        $s =~ s/Cys/C/;
        $s =~ s/Glu/E/;
        $s =~ s/Gln/Q/;
        $s =~ s/Gly/G/;
        $s =~ s/His/H/;
        $s =~ s/Ile/I/;
        $s =~ s/Leu/L/;
        $s =~ s/Lys/K/;
        $s =~ s/Met/M/;
        $s =~ s/Phe/F/;
        $s =~ s/Pro/P/;
        $s =~ s/Ser/S/;
        $s =~ s/Thr/T/;
        $s =~ s/Trp/W/;
        $s =~ s/Tyr/Y/;
        $s =~ s/Val/V/;
    }
    print $s;
}
Ala A
Arg R
Asn N
Asp D
Cys C
Glu E
Gly G
His H
Ile I
Leu L
Lys K
Met M
Phe F
Pro P
Ser S
Thr T
Trp W
Tyr Y
Val V

