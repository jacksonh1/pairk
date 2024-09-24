.. _k-mer-alignment-in-depth:

========================================
PairK alignment - in depth
========================================

there are 2 main ways to run the k-mer alignment step:

-  **Scoring matrix alignment**

   -  These methods use a scoring matrix to score the query k-mer to
      homolog k-mer matches and select the best scoring match from each
      homolog.

-  **Embedding distance alignment**

   -  This method uses the Euclidean distance between the query k-mer
      residue embeddings from a protein large language model (such as
      ESM2) and homolog k-mer residue embeddings and selects the lowest
      distance match from each homolog.

The results from these methods are stored in a :class:`pairk.PairkAln`.
Each of these methods is described in more detail below.

*****************************
Scoring matrix alignment
*****************************

There are 2 implementations of the scoring matrix method:

-  :ref:`pairk.pairk_alignment <1a>` - the original implementation. This is a
   bit slow because it does an exhaustive comparison of all k-mers in
   the query sequence with all k-mers in the homologs. This method is mainly
   included because it is the most flexible and would be the easiest to 
   modify if you wanted to develop a new scoring scheme to compare k-mers.

-  :ref:`pairk.pairk_alignment_needleman <1b>` - a faster implementation that
   uses the Needleman-Wunsch algorithm (as implemented in Biopython) to
   align the k-mers. This is faster and should yield the same results.

inputs common to both functions are:

-  ``idr_dict``: a dictionary of IDR sequences, where the keys are the
   sequence ids and the values are the sequences. Includes the query
   sequence (the sequence to split into k-mers and align with the
   homologs).
-  ``query_id``: a query sequence id (the sequence to split into k-mers
   and align with the homologs). This id should be present in
   ``idr_dict``.
-  ``k``: the length of the k-mers

These inputs can be generated many ways, and there are helper functions
in the pairk library to help with this (see functions in
``pairk.utilities``).

.. _1a:

Scoring matrix - :func:`pairk.pairk_alignment`
================================================================

.. code:: ipython3

    aln_results = pairk.pairk_alignment(
        idr_dict=ex1.idr_dict,
        query_id=ex1.query_id,
        k=5,
    )

To specify the scoring matrix used, you can pass the name of the matrix
to the ``matrix_name`` argument.

To see the available matrices, use the
:func:`pairk.print_available_matrices()` function.

.. code:: ipython3

    pairk.print_available_matrices()


.. parsed-literal::

    biopython-builtin matrices (aligner compatible):
    BENNER22
    BENNER6
    BENNER74
    BLASTN
    BLASTP
    BLOSUM45
    BLOSUM50
    BLOSUM62
    BLOSUM80
    BLOSUM90
    DAYHOFF
    FENG
    GENETIC
    GONNET1992
    HOXD70
    JOHNSON
    JONES
    LEVIN
    MCLACHLAN
    MDM78
    MEGABLAST
    NUC.4.4
    PAM250
    PAM30
    PAM70
    RAO
    RISLER
    SCHNEIDER
    STR
    TRANS
    
    other matrices:
    grantham_similarity_normx100_aligner_compatible
    BLOSUM62
    EDSSMat50
    grantham
    grantham_similarity_norm


.. _1b:

Scoring matrix - :func:`pairk.pairk_alignment_needleman`
================================================================================


This method returns the same results as :func:`pairk.pairk_alignment`, but
it is faster.

The difference is that the :func:`pairk.pairk_alignment_needleman` method
uses the Needleman-Wunsch algorithm to align the k-mers, while the
:func:`pairk.pairk_alignment` method uses a scoring matrix to exhaustively
score the k-mer matches. :func:`pairk.pairk_alignment_needleman` ensures
that the alignment is gapless by using an extremely high gap opening and
extension penalty (-1000000). This will ensure that the alignment is
gapless, unless you use a really unusual scoring matrix with very high
scores.

This method takes similar arguments as :func:`pairk.pairk_alignment`, accept
that the :func:`pairk.pairk_alignment_needleman` method takes an optional
``aligner`` argument. This allows you to create the aligner before
calling the method, which is useful if you want to do multiprocessing,
so that you’re not creating a new aligner for each process. I’ve found
that if you create a new aligner for each process, the memory usage gets
very high.

The ``aligner`` object can be created via the :func:`pairk.create_aligner`
function. This function takes the name of the scoring matrix as an
argument and returns the aligner object. If you don’t pass the
``aligner`` argument to the :func:`pairk.pairk_alignment_needleman` method,
it will create a new aligner using the ``matrix_name`` argument. This is
fine if you’re not doing multiprocessing. If you are doing
multiprocessing, I would suggest creating the aligner before calling the
method. If the ``aligner`` argument is passed, the ``matrix_name``
argument is ignored.

.. code:: ipython3

    aligner = pairk.make_aligner('EDSSMat50')

.. code:: ipython3

    aln_results_needleman = pairk.pairk_alignment_needleman(
        idr_dict=ex1.idr_dict,
        query_id=ex1.query_id,
        k=5,
        aligner=aligner
    )

results are the same as the ``pairk.pairk_alignment`` method in this
case

.. code:: ipython3

    (aln_results.position_matrix == aln_results_needleman.position_matrix).all().all()

.. parsed-literal::

    True


*********************************
Embedding distance alignment
*********************************

This method uses the Euclidean distance between the query k-mer residue
embeddings and homolog k-mer residue embeddings and selects the lowest
distance match from each homolog. For each homolog, it calculates the
distance between the query k-mer and each k-mer in the homolog. It then
selects the k-mer with the lowest distance as the best match.

embedding distance - :func:`pairk.pairk_alignment_embedding_distance`
================================================================================


This method generates residue embeddings using the ESM2 protein large
language model unless residue embeddings are provided directly by the
user.

Because residue embeddings are used, the inputs are slightly different
than the previous methods.

The inputs are:

-  ``full_length_sequence_dict``: a dictionary of full-length sequences,
   where the keys are the sequence ids and the values are the sequences.
   This is used to generate the embeddings.
-  ``idr_position_map``: a dictionary where the keys are the full-length
   sequence ids and the values are the start and end positions of the
   IDR in the full-length sequence (using python indexing). This is used
   to slice out the IDR embeddings/sequences from the full-length
   embeddings/sequences.
-  ``query_id``: a query sequence id (the sequence to split into k-mers
   and align with the homologs). This id should be present in
   ``idr_position_map`` and ``full_length_sequence_dict``.
-  ``k`` - the length of the k-mers
-  ``mod`` - a ``pairk.ESM_Model`` object. This is the ESM2 model used
   to generate the embeddings. The code for the ESM2 embeddings is
   adapted from the kibby conservation tool
   `link <https://github.com/esbgkannan/kibby>`__ DOI:
   10.1093/bib/bbac599
-  ``device`` - whether to use cuda or your cpu for pytorch, should be
   either “cpu” or “cuda”. (default is “cuda”). If “cuda” fails, it will
   default to “cpu”. This argument is passed to the
   ``pairk.ESM_Model.encode`` method

Full length sequences (``full_length_sequence_dict``) are required to
generate the embeddings because each embedding is dependent upon the
neighboring residues. The embeddings for just an IDR are different than
the embeddings for full length sequences. Thus, the full length
embeddings are gathered first, and then the IDR embeddings are sliced
out for the k-mer alignment.

The ``idr_position_map`` is used to slice out the IDR embeddings, and
there must be IDR positions for each sequence in the input sequence set.

Our example has both of those

.. code:: ipython3

    ex1.full_length_dict


.. parsed-literal::

    {'9606_0:00294e': 'MGESSEDIDQMFSTLLGEMDLLTQSLGVDTLPPPDPNPPRAEFNYSVGFKDLNESLNALEDQDLDALMADLVADISEAEQRTIQAQKESLQNQHHSASLQASIFSGAASLGYGTNVAATGISQYEDDLPPPPADPVLDLPLPPPPPEPLSQEEEEAQAKADKIKLALEKLKEAKVKKLVVKVHMNDNSTKSLMVDERQLARDVLDNLFEKTHCDCNVDWCLYEIYPELQIERFFEDHENVVEVLSDWTRDTENKILFLEKEEKYAVFKNPQNFYLDNRGKKESKETNEKMNAKNKESLLEESFCGTSIIVPELEGALYLKEDGKKSWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGTQHKMKYKAPTDYCFVLKHPQIQKESQYIKYLCCDDTRTLNQWVMGIRIAKYGKTLYDNYQRAVAKAGLASRWTNLGTVNAAAPAQPSTGPKTGTTQPNGQIPQATHSVSAVLQEAQRHAETSKDKKPALGNHHDPAVPRAPHAPKSSLPPPPPVRRSSDTSGSPATPLKAKGTGGGGLPAPPDDFLPPPPPPPPLDDPELPPPPPDFMEPPPDFVPPPPPSYAGIAGSELPPPPPPPPAPAPAPVPDSARPPPAVAKRPPVPPKRQENPGHPGGAGGGEQDFMSDLMKALQKKRGNVS',
     '9793_0:005123': 'MGESNEDIDQMFSHLLGEMDLLTKSLGVDTLPPPDPKPPRAEFNFSVGFKDLNESLNALEDQDLDALMADLVADISEAEQRTIQAQRESSQDQLHSASLEKSNFSGAASLGYGADMAIMSTSQYGDELPPPPADPMLDLPLPPPPPEPLSQEEQEAAAKADKIKLALEKLKEAKVKKLVVKVHMNDNSTKSLMVDERQLARDILDNLFEKTHCDCNIDWCLYEIYPELQIERFFEDHENVVEVLSDWTRDTENKVVFLEKEEKYAVFKNPQNFYLDNKGKKESKETNGKMNAKNKESLLEESFCGTSIIVPELEGALYLKEDGKKSWKKRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGTQCKMRYKAPTDYCFVLKHPQIQKESQYIKYLCCDDARTLNQWVTGIRIAKYGKTLYDNYQRAMARAGLASRWTNVGTGNAATPAPPSTGLKTGTAQANGQIPQAAHSVSTVLNEADRQVDTPKDKKPALSNHDPGTPRAQHLPKSSLPPPPPVRRSSDTSSSPVMPAKGAAGGLPPLLDDSLPPPPPPPPLEDDELPPPPPDFDDAPPNFVPPPPPWDAGASLPPPPPPPPPALALAPEATKPSPVVAKRPPVPPKRQENPAPASGGGGGEQDFMSDLMKALQKKRGNVA',
     '1706337_0:000fc7': 'MGDPGACRLDELMRVSTMGESNEDIDQMFSNLLGEMDLLTQSLGVDTVPPPEPKPPRAEFNYSVGFKDLNESLNALEDQDLDALMADLVADIREAEQRTIQAQKESSQSRPHSASLDIPSFSGAASLGYTANVAAPSINQYEDDLPPPPADPMLDLPLPPPTPEPLSQEEEEAQAKADKIKLALEKLKEAKVKKLVVKVHMNDNSTKSLMVDERQLARDVLDNLFEKTHCDCNVDWCLYEIYPELQIERFFEDHENVVEILSDWTRDTENKVLFLEKEEKYAVFKNPQNFYLDNRGKKESKETNEKMNAKNKESLLEESFCGTSIIVPELEGALYLKEDGKKSWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGIQCKMKYKAPTDYCFVLKHPQIQKESQYIKYLCCDDARMLNQWVTGIRIAKYGKTLYDSYQRAMARAGLASRWTNLGTVNAAPPAPSSTGVKTGTTQANGQIPQAAHSMSTVLGEAQRQVETTKDKKSGLGSHDPGAPRAQTLPKSSLPPPPPVRRSSEVGCGSPGTSPKVKGAAAGFPAPPHDLLPPPPPPPPLEDDELPPPPPDFSDAPPDFVPPPPPPSFAGDAGSSLPPPPPPPALAPEAAKPTPVVVKRPPAPPKRQANPGPPGGGGGEQDFMSDLMKALQKKRSNMP',
     '51337_0:001b5a': 'MGESNEDIDQMFSTLLGEMDLLTQSLGVDTLPPPDPKPPRAEFNYSVGFKDLNESLNALEDQDLDALMADLVADISEAEQRTIQAQKESFQNQSHFAPPETSAHSSAACHGDAAHAASITISQCEGDLPPPPADPVLDLPLPPPPPEPLSQEEEEALAKADKIKLALEKLKEAKVKKLVVKVHMFDNSTKSLMVDERQLARDVLDNLFEKTHCDCSVDWCLYETYPELQIERFFEDHENVVEVLSDWTRDTENKVLFLKKEEKYAVFKNPQNFYLDNKGKKENKETSEKMNAKNKEYLLEESFCGTSVLVPELEGALYLKEDGKKSWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGIQCKMKYKAPTDYCFVLKHPQIQKESQYIRHLCCDDAHTLHQWVMGIRIAKYGKTLYDNYQRAVARAGLASRWTNLGTVNTATPAQPSTGFKTGSSQPNGQIPQTIPSVSAGLQEAQRHETIKDKKPSLSSTEPGAPRDPPGARSSLPPPPPPVRRSSDTCARAASPFPAPPDDLPPPPPPPPLEDPAMLPPPPALPEPPPDCVPPPPPPPGPGPQPARPSPGAGRRPPVPPKRQENPGLPSAGAGGEQDFMSDLMKALQKRGHMP',
     '9568_0:004ae1': 'MGESSEDIDQMFSTLLGEMDLLTQSLGVDTLPPPDPNPPRAEFNYSVGFKDLNESLNALEDQDLEALMADLVADISEAEQRTIQAQKESSQNQHHSASLQASNFSGAAPLGHGTNVAATGISQYEDDLPPPPADPVLDLPLPPPPPEPLSQEEEEAQAKADKIKLALEKLKEAKVKKLVVKVHMDDNSTKSLMVDERQLARDVLDNLFEKTHCDCNVDWCLYEIYPELQIEKESVKSVHVRYNLIRGKVSSCKVVPQNFYLDNRGKKESKETNEKMNAKNKESLLEESFCGTSIIVPELEGALYLKEDGKKSWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGTQHKMKYKAPTDYCFVLKHPQIQKESQYIKYLCCDDARTLNQWVMGIRIAKYGKTLYDNYQRAVAKAGLASRWTNLGTVNAAAPAQPSTGPINGTAQPNGQMPQAAHSVSAVLQEAQRHAETSKVKPARPINGTAQPNGQMPQAAHSVSAVLQEAQRHAETSKRPSPAVAKRPPMPPKRHENPGTPSGAGGGEQDFMSDLMKALQKKRGNVS',
     '43346_0:004190': 'MDVQAEVNLPLHKLFFLTYLPIDGPFRGKKEMGQLSPRTRLLLDCLFLDFFLFQMGESNEDIDQMFSNLLGEMDLLTQSLGVDILPPPDPKPPRAEFNYSVGFKDLNESLNALEDQDLDALMADLVADISEAEQRTIQAQKESPQKQSPSASLCVPSFSDTASLGYGANVAAPSQYDDDLPPPPADPMLDLPLPPPPPGPISQEEEEAQAKADKIKLALEKLKEAKVKKLVVKILMNDNSSKSLMVDERQLARDVLDNLFEKTHCDCNVDWCLYEIYPELQIERFFEDHENVVEILSDWTRDTENKLLFLEKEEKYAVFKNPQNFYLDNKGKKENKETNEKMNAKNKESLLEESFCGTSVIVPELEGALYLKEDGKKSWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGIQCKMKYKAPTDYCFVLKHPQIQKESQYIKYLCCDDARALSQWVTGIRIAKYGKTLYDNYQRAMARAGLASRWTNLGTVNAVPPAPPSTGVKTGTTQANGQLPQATQSMNTALGEDWRQLETTKDKKPGPGLGSHDPGAPRAQPLPKSSLPPPPPPVRRSSDVGGAPPPSFAEDLPPPPPPPPALAPESVRTPPVVVKRPPPPPKRQENPGPPGGGGGEQDFMSDLMKALQKKRGNVS',
     '885580_0:00488c': 'MGESHDDIDQMFSTLLGEMDLLTQSLGVDTLPPPAPEPPRAEFNYTVGFKDLNESLNALEDQDLDALMADLVADISEAEQRTIQAQRESCQSQNHASSLGASDCGGVTSVGYAANVTAVGISPYEDDLPPPPDDPMLDLPPPPPPPEPLSKEEEEAQAKADKIKLALEKLKEAKVKKLVVKVHMNDNSTKSLMVDERQLARDVLDNLFEKTHCDCSVDWCLYEIYPELQIERFFEDHENVVEVLSDWTRDTENKVLFLEKEEKYAVFKNPQNFYLDNKGKKEVKQTKEKMNAKTKESLLEESFCGTSVIVPELEGALYLKEDGKKAWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNVYYGVQCRVKYKAPSEHCFVLKHPQIQKESQYIKYLCCDDARTLHQWVTGIRVAKYGKTLYENYQRAVARAGLASRWTNLGTVNTAAAAQPPAGLRTGTSQPNGQLPRDAPCLSDALRESQRQADTSKVGPEGQGGHELPDRPKAAFQASTALAASPYSPKH',
     '10181_0:00305d': 'MGQSHHLFTWTIARESLLSSPLHPSLFSCPMQLALGRITRQSYQMGESNDDIDQMFSTLLGEMDLLTQSLGVDTLPPPDPDPPRPEFNYTVGFKDLNESLNALEDQDLDALMADLVADISEAEQRTIQAQKESSQSQNHAASLEASDCSGDASVGSGANVTAVNISQYEDELPPPPADPMLDLPPPPPPPEPLSKEEEEAQAKADKIRLALEKLKEAKVKKLVVKVHMNDNSTKSLMVDERQLARDVLDNLFEKTHCDCSVDWCLYEIYPELQIERFFEDHENVVEVLSDWTRDTENKVLFLEKEEKYAVFKNPQNFYLDNKGKKEGKKTNEKMNAKNKESLLEESFCGTSIIVPELEGALYLKEDGKKAWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNIYYGVQCKMKYKAPTDHCFVLKHPQIQKESQYIKYLCCDDARTLNQWVTGIRIAKYGKTLHENYQRAVARAGLASRWTNLGTVSTATAAQPPTGLRTSTTQHNGQLPQAAPHLSAVLQEAQRQAEASKVGPEGNRPKLQQISAPSPPSKPEAVFWHLLLFPASENPPYCNFT',
     '1415580_0:000900': 'MEQACDDIDQMFSDLLGEMDLLTQSLGVETIPPPSPKAPNTEFNFSVGFKDLNESLNALEDNDLDALMADLMADINETEKKTFQAQKPPSSSQRSTFTDPEPGFSIAASFDYQSNIPAAYTQDFEDYLPPLPLPPKLDLALPLPPPEPSEPLSKEELESKAKTDKIKLALEKMKEAKVKKLVIKVLMDDDSSKTLMVDERQTARDILDNLFEKTHCDCNIDWCLYEVYQELQIERWFEDHENIVDALSGWTRDSENKMLFLQKKEKYAVFKNPQNFYLAKKGKDAGKEMNEKNKESLLKESFCGTSVIVPELEGALYLKEDGKKSWKKRYFLLRASGIYYVPKGKTKTSRDLACFIQFDNVNIYYGTQYKVKFKAPTDHCFVLKHPQIQKESQYIKYLCCDDSWTLHQWVTGIRIAKFGKTLYDNYRFAVQQMGLASRWPNLSKVDPTVTARSSSSGAVQANGQIHQNVIPVISTNPEAFKRAEDKKPNVGRKPDQAQLQPVPSSNHQQSPKVALHASKIPPPAPARISSQAYSSALTLPSNVKNVNANVLLPPPPSPSPPPPDAFPLPPPCNNDLPPPPDDFYDPPPDFLPPPPPCFATGDRAQLPPGPPLPPPPPSSNQPKPFMKKPVPLPPKRQDITSLHSEQPSLAGPTPVGGGGGQPDFMSDLMKALQKKRGSTS',
     '61221_0:00105a': 'MEEACEDIDQMFSDLLGEMDLLTQSLGVETLPPPPTKASSDEFNFAVGFKDLNESLNALEDTDLDALMADLVADISEVERSTLQEPKDASRYQQVGTMQPSADAGASYGCRTDLSNIANAHVDSGLPPPSVELDFDLPPPPPSTPPAPPSELLTKEEQETQAKADKIKLALEKLKEAKVKKLVIKVHMNDNSTKSLMVDERQVARDVLDNLFEKTHCDCNIDWCLYEMCPELQIERFFEDHENVIGVLSDWTRDSENKLLFLEKSEKYAVFKNPQNFYLSNKGKNEIKVMNEKSKESLLEESFCGTSVIVPEVEGALYLKEDGKKSWKKRYFLLRASGIYYVPKGKTKTSRDLACFIQFENVNVYYGIQYKMKYKAPTDHCFVLKHPQIQRESQYIKYLCCEDQQALHHWVTGIRIAKYGKTLYDNYIRAAHKAGLASRLAKSGNTESIVAVGTSAKGSTHANGQVPQSITLSKTDSSETGKSAEMPKVKKPDNTADSIQAPTPNTPQLKHQRKAGGSHSAPPMPPQRVSSAVTAPLQLPTNAEGKGKVCPSDAAEFPPPPESMLPPPELEDLPLPPPPPPEYFESPPDFIPPPPPSCAVAVSAGAPPLPPPPPSASLPRMPLSIKKKPPPPPRRQEESAGQAGLPKPSAPPPKTETAGQGDFMSDLMKALEKKRGATS',
     '7897_0:0033c5': 'MEQACDDIDAMFSDLLGEMDMLTQSLEVESLHPAVPPALNTDFSFAVGFKDLNESLSALEDTDLDALMADLVADINKVEVETSKGHNSALQDSALPPPPSEDVGVIASSIYIPAFPGSTAVNTGYFSEPLPPPPPLPRPPPDADILLPSKPSETLTQEEIEAKAKADKIKDALEKLKEAKVKKLVIKVHMNDESSKTLMVDERQPVREILDNLFEKTHCDCCVDWCLYEINPDLQIERFFEDHENLVEILLHWTRDSENKILFVEHKEKYAVFRNPQNYYLAKKGKGAEKDMKEKMKESLLEESFCGTSVIVPELEGALYLKEDGKKSWKRRYFLLRASGIYYVPKGKTKTSRDLACFIQFDNVNVYYGMQYKVKYKAPTDHCFILKHPQIQKESQYIKYLCCDNQWTLYQWVTGIRIAKYGKTLYDNYKIAIQKAGLASRWANFAKVQSQNTTASSTGSVQANGHAGQVTPVSVSFSEAWKRGDVGKEKQSNDGQDNLPPPPPPPPPPMQGFMNEAFPPPPSLPPIASGSLPPPLRASASAPAPPPISNNFPPPLDELSPPPDDFDFPEPPPDFLPPPPTVSASGVPPPPPPPPPPPAPTAASQPTPLPKKSVPPRRQENTTLSQPRGGGGGGQPDFMSDLTKALQKKRGNAS',
     '8407_0:002bff': 'MEQTCEEIDEMFSNLLGEMDLLTQSLAVESPPPPTTTKAGTDGNVLFGFKDFDSLNTLEDNDLDALMAELVADCTDAEMKVNNQISNSVAFSSNVDSFAYNIPDFTHDLPTGSTDDLSFLPPPPPSEWEMDLPPPPPDPEEPTEKDALESRDKVDKIMLALDKMKEAKVKKLIVKIHMTDNSTKTLMVDERQTVRDILDNVFEKTHCDCTIEWCLYEENPDLQIERFFEDHENIVEILSDWTRDSENKIRFLKKNEKYAVFKNPQNFYMARKGSADMKDMNEKSKVSLLEQSFCAASVVVPDLEGAIYLKEDGKKSWKKRYFLLRASGIYYVPKGKTKTSRDLACFIQFDNVNVYYGMQYKVRYKAPTDHCFVLKHPQIQKESQYIKYLCCDEPWVLHQWVTGIRIAKYGKVLYENYKSAIHKAGLASRWSSLSSSTTPTGAPQGNGQISQNVANVSSSFSDAWKRGETAKDKQQPTEVRKPEQKISLPPSTKQPPPAPVRRPSNAHVVGTPPLPIKAKPVTSNMPPPPPPAEASQWGDDFLPPPPPPELLDTPPNFLPPPPPSFNSESDYPAPPQFTNVGSAGGPPPPPPPPPPPPAALSPKSAPPQLPVKKLPPKPPMRRDSTGQRPNQQNSLMTNGGGAGGQPDFMSDLMSALQKKRSTTT',
     '173247_0:004550': 'MEDIDAMFSDLLGEMDLLTQSLGQETVPPEALPPTNQEVNYSIGFTDLNESLHELEDNDLDALMADLVADLNATEEKLAAEIKGLKTPSPTTSDLPPPPKGLSLHPPSPHPLPASPASSTSSSVSTPASSAASSLPPPPPQNAKPTKEEIEAQMKADKIKLALEKLKEAKVKKLVVKVLMNDGSSKTLMVDERQNVREVLDNLFEKTHCDCNVDWSLCETNPELQLDRTFEDHENLVEPLSAWMRDSENQVLFQERGEKYEVFKNPQNFYLWKKDKRALKDIKDKDKEILIQENFCGTSIIVPDLEGVLHLKEDGKKSWKQRLFQLRASGIYYVPKGKTKSSRDLVCFLQFDNVNVYYGKDFKTKYKAPTDFCFVLKHPQIQKESQYIKYLCCDDARTMNLWVTGIRIAKYGVSLYENYKTAEKKAAVNSVWMNRSTPSSSNPSTPSPTIKAKTPNQANGHAPKPQPTAPDSMDFGNFPPPPSADILPPPPPDPAFPPPPPSLPAKSSSRPVAPQHKLPANFPPPPMAMDNLPPPPLPPPIDDSPEAPPDFLPPPPPAAGFGSLPPPPLSMNSLPPPPHFGGMDQSLPPPPPDPEFLPPPPPEPVFTGAGAPPPPPPPPPPPPAQAAAVPRAPVRPSGSVRKVPPAPPKRTTPSLQVGGGGGGGDFMSELMVAMQKKRGDH',
     '30732_0:0046dd': 'MEDIDAMFSDLLGEMNLLTQSLGQEAAPAADPPTSTKEVNFSIGFSDLNASLNELEDKDLDDLMSDLMADLNATEEKLAAELQSLEAPPPPDLPPPPKGLIAPAAAAPTSPASPASVCTPSSTATSPLPAPPPQSVKPSKEDLEAQLKAEKIKLALEKLKEAKVKKLVVKVLMNDGSSKTLMVDERQTVREVLDNLFEKTHCDCNVDWSLRETNHELQLERTFEDHENLVEPLSAWTRDSENKVLFLERGEKYEVFKNPQNFYLWKKDKKALKEIKDKDKEILIKENFCGTSTIVPDLESVLHLREDGKKSWKQRLFQLRASGIYYVPKGKTKSSRDLVCFVQFDNINVYYGNDFKTKYKAPTDFCFVLKHPQIQKESQYIKYLCCDDAWTMNLWVTGIRIAKYGSVLYENFKTAEKKAVVSSAWTNRSTPSSSNPSTPSPTVKAKAQSQANGHAPKPQPGPVSQDFGHLPPPPPPCPNDDLPPPPPDPVFPPPPPPLAAKRSPKTAGRSQHPQGNFPPPPPEMDHLPPPPPMEESPPDFLPPPPPMNSLPHPPPPPASFGGVDHSLPPPPPDPEFLPPPPPDPQVTGGGGPPPPPPPPPPPPPASAPAPRGALRPTGSAKKMPPAPPKRTTPVMGGGGGGGGGGGGDFMSELMKAMQKKRSDQ',
     '241271_0:0048e4': 'MDDIDAMFSEMLGEMDLLTQSLDSSLGPETLPPEPLPSTNKEVNYSFGFTDLNASLHELEDNDLDALMADLVADLSATEEKLAAQIEDLKMPSPPPSDLPAPPVGLSTHPTSSIASPTSPASSTSSNVSTPASSSTSPLPPPPPQAAKPTMEEIEAQMKADKIKLALEKLKEAKVKKLVVKVLLNDGSTKTLMVDERQSVREVLDNLFEKTHCDCNVDWSLCETNPELNLERTFEDHENLVEPLLAWTRDSENKILFQERPGKNEVFKNPQNFYLWKKDKRALKEIKDRDKELLVQENFCGTSIIVPDLEAVLYLKEDGKKSWKQRLFQLRASGIYYVPKGKTKSSRDLVCFIQFDNVNVYYGKDFKSKYKAPTDFCFVLKHPQIQKESQYIKYLCCDDAWTMNLWVTGIRIAKYGVSLYENFKAAEKKAANSVWTNRTPVSSNQSTPSPTIKAKSPNQANGHAAKPQPGPESQDFGNIPPPPPPPPPPMTGFLPPPPPDPVLPPPPPLLAAKSPKPSPPQRNLPTNFPPPAMDNLPPPPPPPMDDSFEDPPDFLPPPPPAAGFGSLPPPPPPVNSFPPPPPSAGFGGMGQSLPPPPPDPGFLPPPPPQPMFTGGGTIPPPPPPPPPPTAAPRAPVRPTGSVKKAPPAPPKRTTPSLHGGGGGGGGGGGDFMSELMMAMNKKRGTT',
     '8103_0:0045e4': 'MDDIDAMFSDLLGEMDLLTQSLGQETVPPASLPSTNEEVNLSIGFTDLNESLNELEDTDLDALMADLMADLNATEEKLAAQIEDLKVPSPPPSDLQPPPKGLSIRPVSSLASPTSPASSTGSNVSTSATSPLPPPPPQAAKPTKEEIEAQIKADKIKLALEKLKEAKVKKLVVKVLMNDGSAKTLMVDERQTVREVLDNMFEKTHCDCNVDWSLCETNPELQLERAFEDHENLVEPLLAWTRDSENKVLFQERGEKYEVFKNPQNFYLWKKDMRALKDIKDQDKELLIQENFCGTSIIVPDLEGVMHLKEDGKKSWKPRLFQLRASGIYYVPKGKTKSSRDLVCFVQFDNLNIYYGKDFKGKYKAPTDFCFLLKHPQIQKESQYIKYLCCDDAGSMNLWVTGIRIAKYGASLYDNYKTAEKKAAVSSVWTNRSTPSSSNPSTPSPPIKAKSPGQANGHALKPEPGPVAQDFGHVPPPPPADILPPPADILPPPPPQTFLPPPPPPLAAKSSSKPSLPQRHLPTNFPPPPPAMINLPRPPQPPPTDDASEAPPDFLPPPPPAAGFSPLFPPPPPLNALPLPPPPVSFRVEDRSLPPPPPDPGFLPPPPPMFTGAGAPPPPPPPPPPPRVAVRPAGSVKKRPPAVPKRTTPSLRGGGGGDFMSELALAMNKKRSAH',
     '56723_0:00152f': 'MDDIDAMFTDLLGEMDLLTQSLDQPTVVPEPLPSTTEMNYSIGFTDLNESLHELEDHDLDALMADLVADINATEEKLTAQMKDQKVPPPPSSDLPAPPKGLSTYSASSIASPTSPASSTGSNVSTPASSSASPLPPPPPQSAKPTMEEIEAQMKADKIKLALEKLKEAKVKKLVVKVLLTDGSSKTLMVDERQNVREVLDNLFEKTHCDCNVDWSLCETNYELNLERIFEDHENLVEPLLAWTRDTENKVLFQERTEKNDMFRNPQNFYLWKKDKKALKEIKDRDKELLVQENFCGTSVIVPDLEAVLYLKEDGKKSWKQRLFQLRASGIYYVPKGKTKSSRDLVCFVQFDNMNVYYGKDFKTRYKAPTDFCFVLKHPQIQKDSQYIKYLCCDDAWTMNLWVTGIRIAKYGASLYENFKAAEKAAVSSVWTNRSTPSSSSSSTPSPTIKAKSPSQANGHAPKPQPGPVSQDFGNVPPPPPPMANILPPPRPDAFLPPAPPPLARKNSAKPPPPQRHLPTNFPPPPPAMDNLPPPPPPPPMDDALEAPPDFLPPPPPAAGFGSLPPPPPPSNSFPPPPPPGSFGSMGQSLPPPPPDPGFLPPPPPQPVFTGAGAPPPPPPPPPPPTAAAAPRAPVRPSGSVKKIPPATPKRTTPSLQGGGGGGGGGGGGGGDFMSELMLAMNKKRST',
     '210632_0:004c0c': 'MVDIDAMFSDLLGEMDLLTQSLEQEVAPPKSLPALPSADKEVNFSIGFSDLNESLGELEDNDLDALMADLMADLNATEEKLAAEIEELKVPSLPPANLPAPKNLSIHPASSITSPPSASPAGSSSTLASSSLPPPPPQSVKPTTEEMEAQMKADKIKLALEKLKEAKVKKLVVKVMMNDGSSKTLMVDERQNVREVLDNMFEKTHCDCNVDWSLCETNPELQLERAFEDHENLVDLLSTWMRDSENKVLFQERKEKNEVFKNPQNFYLWKKDKKALKDIKDKDKGLLIQENFCGTSIIVPDLEGILHLKEDGKKSWKQRLFQLRASGIYYVPKGKTKSSRDLVCFVQFDNVNVYYSKDYKSKYKAPTDFCFVLKHPQIQKESQYIKYLCCDDAWTLNLWVTGIRIAKYGVALHENFKTAEQKAATSSAWANRSTSSSSNSSTPSPTIKAKSSSQANGIFPKPGPAPQDFGDLPPPPPLAANILPPPPPEPGLPPPPPPPPPQAAKGSAKPAPPKRQMPANFPPPPTAMDNLPPPPPPPPIDNSEAPPDFLPPPPPASGFGSFPPPPPLNSLPPPPRPGGFGGMDQSLPPPPPDPEFLPPPPPPPQAVFTGGGAPPPPPPPPPPPAAAAPSTAIPRVGLRPAGSLKKLPPAPPKRTTPSMQGSGGGGGGDFMSELMLAMQKKRGDHPPAVLASGT',
     '31033_0:00264e': 'MKMDDIDAMFSEMLGEMDLLTKSLDQEMAPPDAPPSTSEEVSFSIGFPDLNESLQELEDSDLDALMADLVADLNATEQKLAAEIEDLKVPPPPQPHLPPKSRGAVSTSSSCSPSPASSATSPLPVPPPQSVKPSMEEIEAQMKADKIKLALEKLKEAKVKKLVVKVLLNDGSSKTLMVDERQSVRDVLDNLFEKTHCDCNVDWSLCETNAELQLERTFEDHENLVEPLLAWTRDSQNKVFFQERPEKNEVFKNPQNFYLWKKDKKTLQAIKDKDKEILIKENFCGTSIIVPDLEAVLHLREDGKKSWKQRLFQLRASGLYYVPKGKTKSSRDLICFVQFDNLNVYYGKDFRSKYKAPTEFCFVLKHPQIQKESQYIKYLCCDDAWTMNLWVAGIRIAKYGTALHQNYQTALRKAAVTSAWTNCSKPSSDGPPTPPTTIKASPANGHVPKPPPGAAPQDVFPPPPPPMDILPPPPPDPAFPPPPPPLMAKRSPKPSAGHRQAPGDHLPPPPLAPPHDDASEDPPDFLPPPPPSFDSLPPPPPGMSAFPPPPPLLGFSETSQPLPPPPPDPELLLPPPPASMISTGAGAPPPPPPPPPAAAASPRPAPTASGSVRKRPPAPPKRTTPALHGSGGGAGEGAGGGDFMSELMKAMNKKRADHS',
     '63155_0:004c86': 'MDDIDAMFSAMLGEMDLLTQSLGEEKAHPEPHPSSDKQVNFSIGFTDLNESLHELEDNDLDMLMADLMADLNATEEKLAAEIHGLKEPPQPKPDPLPLPRGSSNAPVSEHILPASSGGSGSSNVSTPAPSAACSLPPPPPQCVKPTMDDIEAQAKADKIKLALEKLKEAKVKKLVVKVLMNDGSSKTLMVDERQNVREVLDNMFEKTHCDCNVNWSLCETNPELQLERAFEEHENLVESLSAWIRDSENKVLFQERPEKYEVFKNPQNFYLWKKDKKTLKDIKDKDKELLIQENFCGTSIIVPDLEGVLHLKEDGKKSWKQRLFQLRASGIYYVPKGKTKSSRDLVCFVQFDNMNVYYGKDFKTKHKAPTDFCFVLKHPQIQKDSQYIKYLCCDDAWTMNLWVTGIRIAKYGASLYENYKTAEIKGSNSMWTNRSTPSSSNQSTPSPTVKAKSPNQANGHPPKPQPGPISQAPFPPPPLAEVLPPPPPDPVLPPPPPMPAKSSAKPSPPKRQQQSNFPPPPPELDNLPPPPPPPPTDDTAEAPPDFLPPPPPAVGFGSLPPPPPSFGGVGQSLPPPPPDPQSLPPPPPDPVFIGAGAPPPPPPPPPPPAPGAPVTTLRPAVRPSGSLKKVPPAPPQRNTPSVSGGGGGGGGDFMSELMLAMQKKRGAQ',
     '7994_0:004d71': 'MDDIDAMFTDMLEEMDLLTQSLGAEATEPTPPSKSSSSSSFNSMPEMSNFSIGFTDLNASLNELEDNDLDSLMADLVADLNATEELFAAEKGGVKEPRPPPAVTVPAVHFGSAAPIAPAAPTPSKPKNDVTSCPPAGNTQSLPPPPPASTRPSTDDPEAQKAEKIKLALEKLKEAKVKKLIVKVEITDGSSKTLMVDERQTVRDVMDNLFEKTHCDGNVDWCLCETNPELQTERGFEDHENLVEPLSAWTRDSENKVLFHERKDKYEVFKNPQNFYMWKKDKKSLMDMKDKDKELLLEENFCGTSVIVPDLEGMMYLKEEGKKSWKQRYFLLRASGLYYLPKGKTKSSKDMVCLVQFDNMNVYYCSEYKTKYKAPTDYCFILKHPQIQKESQYIKYLCCDDKWTMTLWVTGIRIAKYGKTMYENYKTAARKGSSLSAVWTSMNRQPSPSTSNTSTPSPTPKAKTANGHASQPRSETVPKAPSNQSAFPPPPPPADFLPPPPPDPTLPPPPPPPPALPVKKESNPPRSAPQRSQPAFPPPPPAMDFSLPPPPPPSDDLEMPPDFLPPPPPAPGGFMGGDLLPPPPPEPFHAPLPPPPAAFHPPPAVHPPPQATGGDLPPPPPPPPPPPPAPAAFHQTPSVRKVGPPPPKRTTPSLAAPSGGDFMSELMLAMNKKRGGQ',
     '109280_0:00369f': 'MDDIDAMFSDLLGEMDLLTQSLGQEQAPPPSSPPEAEQEVNLSIGFTDLNASLNELEDNDLDALMADLVADLNATEEKLAAEIESLKEPQPEPLPPPSVGPPSSSPPLSSDSSTTFPSSTLPPPPPQSSKPTMEEIEAQIKADKIKLALEKLKEAKVKKLVVKVCMNDGSSKTLMVDERQNVREVLDNLCEKTHCDYNVDWSLCETNPELQLERTFEDHEHLVEPLLAWTRDSENQILFQESSDKYEVFKNPQKFYLWKKDKKVLKDMKDKDKEILIKENFCGTSIMVPDLEGVLHLKEDGKKSWKPRLFQLRASGIYYVPKGKTKSSRDLVCFVQFDNVNVYYGKDFRAKHKAPSDFCFVLKHPQIQKDSQYIKYLCCDDAWTMNLWVTGIRIAKYGSNLYENFKTAEKKAAVGSAWASCSVTSGQKQSQSQVANGHANNTSPSPLPPPPPPLGEDLPPPPPPPPQLGKTLPPAPPPLGATLPTPPPPPGGTPPPPPPPRRNPPPYPRHLPHISELYPLRRLLLLYQLPTVHSLVLLFQPNLPPNPSPNHDVQRPISRCLPGPQITFPPPPPPPVDDSPPDFLPPPPPAANFGSHPPPPPPVKTLPPPPPHMKTLPPARLSFKSTNLPPPPPDPGFLPPPLTGVPPPPPPPPPPPPTTAAAGPRRAPVRPSGSLKKMPPPPPKRSTPSLHGRRDGDRGDGDGGGGGDFMSELMRAMQKKRDPH',
     '150288_0:004e5a': 'MDDIDVMFSHLLEEMDDLTQSLVQSADTAADAQTNSSGASDLNEPLNNLDKSESDHLLAQPEETLPVDNQTAPSDPPLPSASSVTLASPTTLAMLKLEPMEVQNESPAPKLTMPQTANNEKPPQTIIKVWMSDGSTKTLMVEGTQTVRDVLDKLFEKTYCDCATEWSLCEINQELHVERILEDHECFVESLSMWSSVTDNKLYFLKRPQKYVIFTQPQFFYMWKRSSLKAISEQEQQLLLKENFGGLTAVVPDLEGWLYLKDDGRKVWKPRYFVLRASGLYYVPKGKTKSSSDLACFVRFEQVNVYSADGHRIRYRAPTDYCFVLKHPCIQKESQYVKFLCCENEDTVLLWVNSIRIAKYGTVLYENYKTALKRAQHPPDRCSTSSDNLNSQIGQSAPTPDECIEQDEPPPDFIPPPPPGYMAIL'}


.. code:: ipython3

    ex1.idr_position_map


.. parsed-literal::

    {'9606_0:00294e': [440, 665],
     '9793_0:005123': [440, 657],
     '1706337_0:000fc7': [457, 676],
     '51337_0:001b5a': [440, 632],
     '9568_0:004ae1': [426, 564],
     '43346_0:004190': [492, 656],
     '885580_0:00488c': [440, 524],
     '10181_0:00305d': [484, 578],
     '1415580_0:000900': [439, 679],
     '61221_0:00105a': [440, 678],
     '7897_0:0033c5': [442, 653],
     '8407_0:002bff': [430, 663],
     '173247_0:004550': [432, 680],
     '30732_0:0046dd': [425, 663],
     '241271_0:0048e4': [435, 685],
     '8103_0:0045e4': [429, 673],
     '56723_0:00152f': [430, 685],
     '210632_0:004c0c': [429, 685],
     '31033_0:00264e': [420, 658],
     '63155_0:004c86': [431, 667],
     '7994_0:004d71': [439, 676],
     '109280_0:00369f': [418, 723],
     '150288_0:004e5a': [376, 424]}


The ``mod`` input is required so that you can preload the ESM model
before running the method. You preload the ESM model with
``pairk.ESM_Model()``


you can use pre-generated embeddings by providing them in a dictionary
format to the ``precomputed_embeddings`` argument. The keys should be
the sequence ids and the values should be full length sequence embedding
tensors. For each sequence, the tensor.shape[0] should be equal to
``l``\ +2, where ``l`` is the length of the full length sequence. The +2
is for the start and end tokens. If you provide precomputed embeddings,
the ``mod`` and ``device`` arguments are ignored.

The ``pairk.pairk_alignment_embedding_distance`` method returns a
``PairkAln`` object, just like the previous methods

example usage: loading the ESM2 model and running the method

.. code:: ipython3

    mod = pairk.ESM_Model(threads=4)
    aln_results_embedding = pairk.pairk_alignment_embedding_distance(
        full_length_sequence_dict=ex1.full_length_dict,
        idr_position_map=ex1.idr_position_map,
        query_id=ex1.query_id,
        k=5,
        mod=mod,
        device="cpu"
    )

****************************
k-mer alignment results
****************************

The results of the above pairwise k-mer alignment methods are returned
as a :class:`pairk.PairkAln` object.


The actual “alignments” are stored as matrices in the :class:`pairk.PairkAln`
object. The main matrices are:

-  orthokmer_matrix - the best matching k-mers from each homolog for
   each query k-mer
-  position_matrix - the positions of the best matching k-mers in the
   homologs
-  score_matrix - the scores of the best matching k-mers

Each matrix is a pandas DataFrame where the index is the start position
of the k-mer in the query sequence. The columns are the query k-mers +
the homolog sequence ids.

The :class:`pairk.PairkAln` object has some useful methods for accessing the data.
For example, you can get the best matching k-mers for a query k-mer by
its position in the query sequence using the :func:`pairk.get_pseudo_alignment`
method (or by directly accessing the dataframes). You can also plot the
matrices as heatmaps, save the results to a json file, and load the
results from that file

example: accessing the DataFrames from the :class:`pairk.PairkAln` object directly

.. code:: ipython3

    aln_results.score_matrix




.. raw:: html

    <div style="max-width: 100%; overflow-x: auto;">
      <style scoped>
          .dataframe tbody tr th:only-of-type {
              vertical-align: middle;
          }
      
          .dataframe tbody tr th {
              vertical-align: top;
          }
      
          .dataframe thead th {
              text-align: right;
          }
      </style>
      <table border="1" class="dataframe">
        <thead>
          <tr style="text-align: right;">
            <th></th>
            <th>query_kmer</th>
            <th>9793_0:005123</th>
            <th>1706337_0:000fc7</th>
            <th>51337_0:001b5a</th>
            <th>9568_0:004ae1</th>
            <th>43346_0:004190</th>
            <th>885580_0:00488c</th>
            <th>10181_0:00305d</th>
            <th>1415580_0:000900</th>
            <th>61221_0:00105a</th>
            <th>...</th>
            <th>30732_0:0046dd</th>
            <th>241271_0:0048e4</th>
            <th>8103_0:0045e4</th>
            <th>56723_0:00152f</th>
            <th>210632_0:004c0c</th>
            <th>31033_0:00264e</th>
            <th>63155_0:004c86</th>
            <th>7994_0:004d71</th>
            <th>109280_0:00369f</th>
            <th>150288_0:004e5a</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th>0</th>
            <td>TNLGT</td>
            <td>22.0</td>
            <td>28.0</td>
            <td>28.0</td>
            <td>28.0</td>
            <td>28.0</td>
            <td>28.0</td>
            <td>28.0</td>
            <td>10.0</td>
            <td>9.0</td>
            <td>...</td>
            <td>13.0</td>
            <td>12.0</td>
            <td>13.0</td>
            <td>13.0</td>
            <td>9.0</td>
            <td>8.0</td>
            <td>13.0</td>
            <td>8.0</td>
            <td>14.0</td>
            <td>9.0</td>
          </tr>
          <tr>
            <th>1</th>
            <td>NLGTV</td>
            <td>16.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>16.0</td>
            <td>8.0</td>
            <td>...</td>
            <td>11.0</td>
            <td>13.0</td>
            <td>11.0</td>
            <td>13.0</td>
            <td>7.0</td>
            <td>7.0</td>
            <td>11.0</td>
            <td>6.0</td>
            <td>14.0</td>
            <td>10.0</td>
          </tr>
          <tr>
            <th>2</th>
            <td>LGTVN</td>
            <td>16.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>23.0</td>
            <td>10.0</td>
            <td>7.0</td>
            <td>...</td>
            <td>11.0</td>
            <td>8.0</td>
            <td>9.0</td>
            <td>10.0</td>
            <td>7.0</td>
            <td>8.0</td>
            <td>10.0</td>
            <td>8.0</td>
            <td>16.0</td>
            <td>4.0</td>
          </tr>
          <tr>
            <th>3</th>
            <td>GTVNA</td>
            <td>19.0</td>
            <td>26.0</td>
            <td>23.0</td>
            <td>26.0</td>
            <td>26.0</td>
            <td>23.0</td>
            <td>17.0</td>
            <td>15.0</td>
            <td>10.0</td>
            <td>...</td>
            <td>11.0</td>
            <td>9.0</td>
            <td>10.0</td>
            <td>8.0</td>
            <td>8.0</td>
            <td>9.0</td>
            <td>11.0</td>
            <td>9.0</td>
            <td>9.0</td>
            <td>6.0</td>
          </tr>
          <tr>
            <th>4</th>
            <td>TVNAA</td>
            <td>18.0</td>
            <td>25.0</td>
            <td>22.0</td>
            <td>25.0</td>
            <td>22.0</td>
            <td>22.0</td>
            <td>16.0</td>
            <td>14.0</td>
            <td>11.0</td>
            <td>...</td>
            <td>11.0</td>
            <td>9.0</td>
            <td>8.0</td>
            <td>12.0</td>
            <td>8.0</td>
            <td>11.0</td>
            <td>11.0</td>
            <td>10.0</td>
            <td>11.0</td>
            <td>5.0</td>
          </tr>
          <tr>
            <th>...</th>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
            <td>...</td>
          </tr>
          <tr>
            <th>217</th>
            <td>LQKKR</td>
            <td>31.0</td>
            <td>31.0</td>
            <td>19.0</td>
            <td>31.0</td>
            <td>31.0</td>
            <td>3.0</td>
            <td>9.0</td>
            <td>31.0</td>
            <td>24.0</td>
            <td>...</td>
            <td>26.0</td>
            <td>19.0</td>
            <td>19.0</td>
            <td>19.0</td>
            <td>26.0</td>
            <td>19.0</td>
            <td>26.0</td>
            <td>19.0</td>
            <td>26.0</td>
            <td>0.0</td>
          </tr>
          <tr>
            <th>218</th>
            <td>QKKRG</td>
            <td>29.0</td>
            <td>23.0</td>
            <td>16.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>5.0</td>
            <td>3.0</td>
            <td>29.0</td>
            <td>22.0</td>
            <td>...</td>
            <td>23.0</td>
            <td>22.0</td>
            <td>16.0</td>
            <td>16.0</td>
            <td>29.0</td>
            <td>17.0</td>
            <td>29.0</td>
            <td>22.0</td>
            <td>23.0</td>
            <td>-1.0</td>
          </tr>
          <tr>
            <th>219</th>
            <td>KKRGN</td>
            <td>29.0</td>
            <td>23.0</td>
            <td>18.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>2.0</td>
            <td>2.0</td>
            <td>23.0</td>
            <td>21.0</td>
            <td>...</td>
            <td>17.0</td>
            <td>23.0</td>
            <td>15.0</td>
            <td>17.0</td>
            <td>23.0</td>
            <td>18.0</td>
            <td>21.0</td>
            <td>22.0</td>
            <td>14.0</td>
            <td>-1.0</td>
          </tr>
          <tr>
            <th>220</th>
            <td>KRGNV</td>
            <td>29.0</td>
            <td>19.0</td>
            <td>20.0</td>
            <td>29.0</td>
            <td>29.0</td>
            <td>8.0</td>
            <td>8.0</td>
            <td>17.0</td>
            <td>15.0</td>
            <td>...</td>
            <td>9.0</td>
            <td>17.0</td>
            <td>8.0</td>
            <td>11.0</td>
            <td>14.0</td>
            <td>9.0</td>
            <td>13.0</td>
            <td>14.0</td>
            <td>11.0</td>
            <td>2.0</td>
          </tr>
          <tr>
            <th>221</th>
            <td>RGNVS</td>
            <td>23.0</td>
            <td>12.0</td>
            <td>13.0</td>
            <td>27.0</td>
            <td>27.0</td>
            <td>9.0</td>
            <td>13.0</td>
            <td>15.0</td>
            <td>13.0</td>
            <td>...</td>
            <td>10.0</td>
            <td>11.0</td>
            <td>8.0</td>
            <td>13.0</td>
            <td>7.0</td>
            <td>9.0</td>
            <td>7.0</td>
            <td>9.0</td>
            <td>8.0</td>
            <td>9.0</td>
          </tr>
        </tbody>
      </table>
      <p>222 rows × 23 columns</p>
    </div>


example: access the best matching k-mers for the query k-mer at position
4:

.. code:: ipython3

    print(aln_results.orthokmer_matrix.loc[4])


.. parsed-literal::

    query_kmer          TVNAA
    9793_0:005123       TGNAA
    1706337_0:000fc7    TVNAA
    51337_0:001b5a      TVNTA
    9568_0:004ae1       TVNAA
    43346_0:004190      TVNAV
    885580_0:00488c     TVNTA
    10181_0:00305d      TVSTA
    1415580_0:000900    NVNAN
    61221_0:00105a      AVSAG
    7897_0:0033c5       TVSAS
    8407_0:002bff       SQNVA
    173247_0:004550     TPNQA
    30732_0:0046dd      TVKAK
    241271_0:0048e4     PVNSF
    8103_0:0045e4       PLNAL
    56723_0:00152f      TAAAA
    210632_0:004c0c     TIKAK
    31033_0:00264e      TIKAS
    63155_0:004c86      TVKAK
    7994_0:004d71       TSNTS
    109280_0:00369f     TTAAA
    150288_0:004e5a     NLNSQ
    Name: 4, dtype: object


example: access the best matching k-mers for the query k-mer at position
4 using the ``get_pseudo_alignment`` method. (the returned list includes
the query k-mer sequence)

.. code:: ipython3

    aln_results.get_pseudo_alignment(4)


.. parsed-literal::

    ['TVNAA',
     'TGNAA',
     'TVNAA',
     'TVNTA',
     'TVNAA',
     'TVNAV',
     'TVNTA',
     'TVSTA',
     'NVNAN',
     'AVSAG',
     'TVSAS',
     'SQNVA',
     'TPNQA',
     'TVKAK',
     'PVNSF',
     'PLNAL',
     'TAAAA',
     'TIKAK',
     'TIKAS',
     'TVKAK',
     'TSNTS',
     'TTAAA',
     'NLNSQ']


you can search for a specific kmer to get its positions. You can then
use the positions to query the matrices.

.. code:: ipython3

    aln_results.find_query_kmer_positions('LPPPP')




.. parsed-literal::

    [75, 113, 127, 157]


.. code:: ipython3

    aln_results.get_pseudo_alignment(75)


.. parsed-literal::

    ['LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'PPMPP',
     'LPPPP',
     'LPDRP',
     'APSPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'LPPPP',
     'IPPPP']



.. code:: ipython3

    aln_results.orthokmer_matrix.loc[[75, 113, 127, 157]].T




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>75</th>
          <th>113</th>
          <th>127</th>
          <th>157</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>query_kmer</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>9793_0:005123</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>1706337_0:000fc7</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>51337_0:001b5a</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>9568_0:004ae1</th>
          <td>PPMPP</td>
          <td>PPMPP</td>
          <td>PPMPP</td>
          <td>PPMPP</td>
        </tr>
        <tr>
          <th>43346_0:004190</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>885580_0:00488c</th>
          <td>LPDRP</td>
          <td>LPDRP</td>
          <td>LPDRP</td>
          <td>LPDRP</td>
        </tr>
        <tr>
          <th>10181_0:00305d</th>
          <td>APSPP</td>
          <td>APSPP</td>
          <td>APSPP</td>
          <td>APSPP</td>
        </tr>
        <tr>
          <th>1415580_0:000900</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>61221_0:00105a</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>7897_0:0033c5</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>8407_0:002bff</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>173247_0:004550</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>30732_0:0046dd</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>241271_0:0048e4</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>8103_0:0045e4</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>56723_0:00152f</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>210632_0:004c0c</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>31033_0:00264e</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>63155_0:004c86</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>7994_0:004d71</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>109280_0:00369f</th>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
          <td>LPPPP</td>
        </tr>
        <tr>
          <th>150288_0:004e5a</th>
          <td>IPPPP</td>
          <td>IPPPP</td>
          <td>IPPPP</td>
          <td>IPPPP</td>
        </tr>
      </tbody>
    </table>
    </div>

|

.. Note:: the k-mers are defined by position rather than sequence. You
    could easily make a variant of this method that uses the unique
    sequences instead. It would make the method slightly faster. The reason
    that I didn’t do this is because I wanted to mimic the LLM embedding
    version of Pairk, where identical k-mers have different embeddings and
    thus are treated as different k-mers.Inclusion of duplicate k-mers does
    alter the final z-scores, so it’s something to be aware of.

example: plot a heatmap of the matrices

.. code:: ipython3

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3,3))
    aln_results.plot_position_heatmap(ax)
    ax.xaxis.set_visible(False)

.. image:: notebooks_rst/pairk_tutorial_simplified_files/pairk_tutorial_simplified_62_0.png


example: save the results to a file using ``write_to_file`` and load
them back into python using ``from_file``:

.. code:: ipython3

    aln_results.write_to_file('./aln_results.json')
    aln_results = pairk.PairkAln.from_file('./aln_results.json')
    print(aln_results)


.. parsed-literal::

    PairkAln object for 222 query k-mers
    query sequence: TNLGTVNAAAPAQPSTGPKTGTTQPNGQIPQATHSVSAVLQEAQRHAETSKDKKPALGNHHDPAVPRAPHAPKSSLPPPPPVRRSSDTSGSPATPLKAKGTGGGGLPAPPDDFLPPPPPPPPLDDPELPPPPPDFMEPPPDFVPPPPPSYAGIAGSELPPPPPPPPAPAPAPVPDSARPPPAVAKRPPVPPKRQENPGHPGGAGGGEQDFMSDLMKALQKKRGNVS
    k-mer length: 5
    

