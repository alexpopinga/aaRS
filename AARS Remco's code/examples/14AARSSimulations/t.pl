$s = 'VGEALDRPSTIQNCYMWHKF';
for ($i = 0; $i < length($s); $i++) {
	$c = substr($s, $i, 1);
	for ($j = 1; $j < 8; $j++) {
		print "AARS_".$c."_$j=$i,";
	}
	print "\n";
}
exit;
print '             <distribution id="'.$c.'.prior" monophyletic="true" spec="beast.math.distributions.MRCAPrior" tree="@Tree.t:AARS_14_1">
                <taxonset id="'.$c.'" spec="TaxonSet">
		            <taxon id="AARS_'.$c.'_1"/>
		            <taxon id="AARS_'.$c.'_2"/>
		            <taxon id="AARS_'.$c.'_3"/>
		            <taxon id="AARS_'.$c.'_4"/>
		            <taxon id="AARS_'.$c.'_5"/>
		            <taxon id="AARS_'.$c.'_6"/>
		            <taxon id="AARS_'.$c.'_7"/>
                </taxonset>
            </distribution>
';

	
