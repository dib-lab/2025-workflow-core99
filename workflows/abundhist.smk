ISECT_METAG = {
    's__Lactobacillus amylovorus': ['SRR8960379',
                                    'SRR8960968',
                                    'SRR12795775',
                                    'SRR10209683',
                                    'SRR17241521',
                                    'ERR3211783',
                                    'ERR8314733',
                                    'ERR8314788',
                                    'ERR1135238',
                                    'SRR17241595',
                                    'SRR11124923',
                                    'SRR8960492',
                                    'SRR11125534',
                                    'SRR17241485',
                                    'SRR14812372'],
    's__Mogibacterium_A kristiansenii': ['SRR14812372',
                                         'SRR5240744',
                                         'SRR17241663',
                                         'SRR8960625',
                                         'SRR12795793',
                                         'SRR11183784',
                                         'SRR17241609',
                                         'SRR8960444',
                                         'SRR11125751',
                                         'SRR8960379'],
    's__Prevotella sp002251295': ['SRR8960491',
                                   'SRR5241534',
                                   'ERR1135213',
                                   'ERR3211839',
                                   'ERR3211879',
                                   'ERR3211876',
                                   'SRR14369134',
                                   'ERR3211783',
                                   'SRR14369306',
                                   'SRR11489781'],
    's__Prevotella sp000434975': ['SRR8960379',
                                   'SRR8960101',
                                   'SRR14812375',
                                   'ERR3211876',
                                   'SRR14812372',
                                   'ERR1135303',
                                   'SRR11489769',
                                   'ERR3211839',
                                   'ERR8314767',
                                   'ERR1135349'],
    's__Cryptobacteroides sp900546925': ['SRR17241596',
                                          'ERR3211812',
                                          'ERR3211974',
                                          'ERR8314767',
                                          'ERR8314788',
                                          'ERR3211971',
                                          'SRR10209683',
                                          'SRR11551383',
                                          'SRR14369219',
                                          'SRR17241485'],
    's__Sodaliphilus sp004557565': ['SRR11125534',
                                     'SRR11805642',
                                     'ERR8314743',
                                     'ERR1135349',
                                     'ERR3211891',
                                     'ERR8314788',
                                     'ERR3211929',
                                     'ERR1135303',
                                     'SRR10209683',
                                     'SRR16235699'],
    's__Colivicinus sp002299675': ['SRR5241534',
                                   'ERR3211906',
                                   'SRR5240744',
                                   'ERR1135346',
                                   'SRR17241663',
                                   'SRR10209667',
                                   'SRR12795793',
                                   'ERR8314733',
                                   'SRR11489783',
                                   'SRR11551383'],
    's__Ornithospirochaeta sp022785155': ['ERR1135273',
                                          'ERR3211801',
                                          'SRR18048904',
                                          'ERR8314788',
                                          'SRR8960625',
                                          'ERR1135213',
                                          'ERR3211929',
                                          'ERR1135304',
                                          'SRR11489781',
                                          'SRR11489783']
}

ISECT_FILES = []
for s, mm in ISECT_METAG.items():
    for m in mm:
#        ISECT_FILES.append(f'outputs.abundhist/{s}.x.{m}.abundhist.png')
        ISECT_FILES.append(f'outputs.abundhist/{s}.x.{m}.isect.sig.zip')

rule do_abundhist:
    input:
        ISECT_FILES

rule do_update_mapping_parquet:
    shell: """
        scripts/csv-to-parquet.py outputs.cds/cds3-genes/mapping-coverage.csv -o inputs.abundhist/ --redo
    """
    

rule abundhist_isect_wc:
    input:
        sig='inputs.abundhist/{s}.cds3.sig.zip',
        metag='/group/ctbrowngrp5/wort/wort-sra/sigs/{m}.sig',
    output:
        'outputs.abundhist/{s}.x.{m}.abundhist.png'
    shell: """
        sourmash scripts abundhist --figure {output:q} -I {input.sig:q} \
            {input.metag:q} -k 21 --bins 50 --ymax=100 \
            --figure-title {wildcards.s:q}.x.{wildcards.m:q}
    """

rule sig_isect_wc:
    input:
        sig='inputs.abundhist/{s}.cds3.sig.zip',
        metag='/group/ctbrowngrp5/wort/wort-sra/sigs/{m}.sig',
    output:
        'outputs.abundhist/{s}.x.{m}.isect.sig.zip',
    shell: """
        sourmash sig intersect -k 21 {input.metag:q} {input.sig:q} \
            -o {output:q} -A {input.metag:q}
    """
