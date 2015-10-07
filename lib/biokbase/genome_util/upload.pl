use strict;
use Data::Dumper;
use JSON;
use Bio::KBase::workspace::Client;
use Bio::KBase::AuthToken;
my $token  = $ENV{KB_AUTH_TOKEN};
if (!$token){
  my $to = Bio::KBase::AuthToken->new();
  $token = $to->{token};
}

#usage perl upload.pl blast.json objname wsname ws_url

my $file = $ARGV[0];
my $id = $ARGV[1];
my $ws_id=$ARGV[2];
my $ws_url = $ARGV[3];



open (FILE, $file) or die ("nnc");
my $ws_doc1 = from_json(join ("", <FILE>));
close (FILE);

my $ws_doc=$ws_doc1;
my $wsc = Bio::KBase::workspace::Client->new($ws_url, token=>$token );
my $metadata = $wsc->save_object({id =>$id, type =>"GenomeUtil.BlastOutput", data => $ws_doc, workspace => $ws_id});
