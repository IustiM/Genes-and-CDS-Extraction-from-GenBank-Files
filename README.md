# Genes-and-CDS-Extraction-from-GenBank-Files
Ghid pentru Utilizarea Scriptului Perl pentru Extragerea Genelor si CDS-urilor

 I. Instalare Modul Perl

Scriptul foloseste modulul Bio::SeqIO din BioPerl. Pentru a instala acest modul deschide CMD si executa urmatoarea comanda pentru a instala BioPerl:

     cpan Bio::Perl
     

 II. Pregatirea Fisierelor
1. Pregateste Fisierul de Intrare:
Pune fisierul GenBank (Homo_sapiens.GRCh38.111.chromosome.19.dat) in acelasi director cu scriptul Perl sau specifica calea completa catre fisier.

 III. Utilizarea Scriptului

   In CMD, ruleaza scriptul cu comanda:

     perl “Nume_script”.pl "Homo_sapiens.GRCh38.111.chromosome.19.dat" "output.txt" 100000 320000

     Argumentele din comanda sunt:
     - “Nume_script”.pl: Numele fisierului tau
     - Homo_sapiens.GRCh38.111.chromosome.19.dat: Fisierul de intrare.
     - output.txt: Fisierul de iesire in care vor fi scrise genele si CDS-urile.
     - 100000: Pozitia de start (numerica).
     - 320000: Pozitia de end (numerica).


IV. Pasul Final

1. Salveaza scriptul si fisierul de intrare in acelasi director.
2. Urmeaza pasii de mai sus pentru a rula scriptul din CMD.
3. Verifica fisierul output.txt pentru a vedea rezultatele.

Acum ai toate informatiile necesare pentru a folosi acest script Perl pentru a extrage gene si CDS din fisiere GenBank.
