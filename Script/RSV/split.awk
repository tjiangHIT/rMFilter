#printf(" %s %s\n", $1,$2 ) > "today_rpt1"
BEGIN {
    OFS="\t"
}
{
    printf("%s\n", $0) >  $1
}
