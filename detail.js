// Detail rendering + locus graphic. Depends on globals from index.html.
const UCSC={hs1:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&position=",
  hg38:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=",
            t2t:"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&position="};
let GENES=null;
async function geneModels(){
  if(GENES===null){
    try{const r=await fetch("data/gene_models.json.gz");
      GENES=r.ok?JSON.parse(await new Response(r.body.pipeThrough(new DecompressionStream("gzip"))).text()):{};}
    catch(e){GENES={};}
  }
  return GENES;
}
function coordOf(d,asm){return (d.coord||[]).find(c=>c.assembly===asm);}

function render(d,h){
  if(!d){$("view").innerHTML='<div class="empty">locus not found in bundle</div>';return;}
  const hg=coordOf(d,"hg38"),t2=coordOf(d,"t2t"),g=d.group_info||{},st=d.structure||{},xg=d.crossgenome||{};
  // graphic assembly: hg38 when available, else t2t. The t2t fallback is what gives the
  // ~2,945 loci absent from hg38 a graphic at all; their lanes are repeats-only.
  const gco=hg||t2, gasm=hg?"hg38":(t2?"t2t":null);
  const links=[hg?'<a class="ucsc" target="_blank" href="'+UCSC.hg38+hg.chrom+":"+(hg.start+1-1000)+"-"+(hg.end+1000)+'">UCSC hg38 ±1 kb</a>':"",
               t2?'<a class="ucsc" target="_blank" href="'+UCSC.t2t+t2.chrom+":"+(t2.start+1-1000)+"-"+(t2.end+1000)+'">UCSC T2T (hs1) ±1 kb</a>':""].join("");
  const dfb=(d.dfam_best||[])[0]||{};
  $("view").innerHTML=
   '<div class="panel"><div class="idline"><span class="cid">'+esc(d.combined_id)+'</span>'+
     '<span class="uid">'+esc(d.uid)+'</span><span class="uid">'+esc(d.versioned_id)+'</span></div>'+
     '<div class="note">Cite the versioned_id. locus_uid is the immutable primary key; combined_id is positional and may re-letter.</div>'+
     '<div style="margin-top:9px">'+links+'</div></div>'+
   '<div class="two">'+
     '<div class="panel"><h2>Locus</h2><dl class="kv">'+
       kv("clade",d.group)+kv("band",d.band)+kv("origin",d.origin)+
       kv("structure",st.structure)+kv("category",st.category)+
       kv("LTR names",st.ltr_names)+kv("internal names",st.int_names)+
       kv("tandem/nested",st.is_tandem_or_nested?"yes":"no")+
       kv("segments",d.segments.length+(d.segments.length?"":" (none stored)"))+
     '</dl></div>'+
     '<div class="panel"><h2>Clade — '+esc(d.group)+'</h2><dl class="kv">'+
       kv("superfamily",g.superfamily)+kv("HERV class",g.herv_class)+
       kv("loci in clade",g.n_loci)+kv("int model",g.intModel)+
       kv("RepBase class",g.repbase_class)+kv("HERVd family",g.hervd_family)+
       kv("dominant LTR",g.dominant_ltr)+
       kv("with flanking LTR",g.frac_with_flanking_ltr==null?null:(100*g.frac_with_flanking_ltr).toFixed(1)+"%")+
       kv("extension verdict",g.extension_verdict)+
     '</dl></div></div>'+
   '<div class="panel"><h2>Locus map ±1 kb ('+(gasm||"—")+')'+
     (gco&&gasm!=="hg38"?' <span class="note" style="font-weight:400">— hg38 coordinate absent;'+
       ' gene models, gEVE ORFs and HERVarium domains are hg38-only and are omitted</span>':"")+
     '</h2>'+
     (gco?'<label class="note" style="display:block;margin:-4px 0 6px">'+
       '<input type="checkbox" id="repall"> show all repeat classes'+
       ' (simple repeats, low complexity, satellites)</label>':"")+
     '<div id="gfx">'+
     (gco?'<div class="note">rendering…</div>'
        :'<div class="note">no hg38 or t2t coordinate — graphic unavailable</div>')+'</div>'+
     '<div class="lg"><span><i style="background:var(--ltr)"></i>LTR segment</span>'+
     '<span><i style="background:var(--int)"></i>internal segment</span>'+
     '<span><i style="background:var(--orf)"></i>gEVE ORF</span>'+
     '<span><i style="background:var(--dom)"></i>HERVarium domain</span>'+
     '<span><i style="background:#b5651d"></i>LINE</span>'+
     '<span><i style="background:#2e8b8b"></i>SINE/Alu</span>'+
     '<span><i style="background:#3b6ea5"></i>LTR (RepeatMasker)</span>'+
     '<span><i style="background:#6b7f3a"></i>DNA</span>'+
     '<span><i style="background:var(--gene)"></i>gene exon (thick) / intron (thin)</span></div></div>'+
   '<div class="panel"><h2>Coordinates</h2>'+coordTable(d,xg)+'</div>'+
   '<div class="panel"><h2>Aliases — '+d.aliases.length+' rows</h2>'+aliasTable(d.aliases)+'</div>'+
   '<div class="panel"><h2>Dfam best alignment</h2>'+
     (dfb.consensus_name?'<dl class="kv">'+kv("consensus",dfb.consensus_name)+kv("accession",dfb.dfam_accession)+
       kv("% identity",dfb.pct_identity==null?null:Number(dfb.pct_identity).toFixed(1))+
       kv("consensus coverage",dfb.cons_cov==null?null:(100*dfb.cons_cov).toFixed(1)+"%")+
       kv("SW score",dfb.sw_score)+kv("quality",dfb.aln_quality)+'</dl>'
      :'<div class="note">not screened against Dfam consensus</div>')+'</div>'+
   '<div class="panel"><h2>gEVE ORFs — '+d.geve.length+'</h2>'+tbl(d.geve,
      ["geve_orf_id","orf_class","hg38_chrom","hg38_orf_start","hg38_orf_end","orf_strand"])+'</div>'+
   '<div class="panel"><h2>HERVarium domains — '+d.domains.length+'</h2>'+tbl(d.domains,
      ["gene","domain_desc","element","status","domain_score","hg38_start","hg38_end","strand"])+'</div>'+
   '<div class="panel"><h2>Overlapping genes — '+d.genes.length+' rows</h2>'+tbl(d.genes,
      ["genome","source","ref_gene_name","ref_gene_id","gene_type","overlap_type","overlap_bp","exon_overlap_bp","n_transcripts"])+'</div>'+
   '<div class="panel"><h2>Segments — '+d.segments.length+'</h2>'+tbl(d.segments,
      ["seg_index","segment_class","repName","repFamily","repClass","chrom","start","end","strand","span","rmsk_sw_score"])+'</div>';
  if(gco) drawLocus(d,gco,gasm);
  const rc=$("repall"); if(rc) rc.onchange=()=>{ if(gco) drawLocus(d,gco,gasm); };
}
function kv(k,v){return "<dt>"+esc(k)+"</dt><dd>"+fmt(v)+"</dd>";}
function tbl(rows,cols){
  if(!rows||!rows.length)return '<div class="note">none</div>';
  const use=cols.filter(c=>rows.some(r=>r[c]!=null&&r[c]!==""));
  return "<table><tr>"+use.map(c=>"<th>"+esc(c)+"</th>").join("")+"</tr>"+
    rows.map(r=>"<tr>"+use.map(c=>'<td class="'+(typeof r[c]==="number"?"mono":"")+'">'+fmt(r[c])+"</td>").join("")+"</tr>").join("")+"</table>";
}
function coordTable(d,xg){
  const rows=(d.coord||[]).map(c=>{
    const u=UCSC[c.assembly]; const pos=c.chrom+":"+(c.start+1)+"-"+c.end;
    return "<tr><td>"+esc(c.assembly)+'</td><td class=mono>'+
      (u?'<a target="_blank" href="'+u+pos+'">'+esc(pos)+"</a>":esc(pos))+
      "</td><td>"+esc(c.strand)+"</td><td class=mono>"+bp(c.span)+"</td><td>"+esc(c.contig_type||"—")+
      "</td><td>"+(c.coord_unreliable?"unreliable":"ok")+"</td><td>"+esc(c.lift_method||"—")+"</td></tr>";}).join("");
  return "<table><tr><th>assembly</th><th>position (1-based display)</th><th>strand</th><th>span</th>"+
    "<th>contig</th><th>flag</th><th>lift</th></tr>"+rows+"</table>"+
    (xg.pct_identity!=null?'<div class="note">cross-genome: '+Number(xg.pct_identity).toFixed(2)+
      "% identity, "+esc(xg.xg_class||"")+", edit distance "+fmt(xg.edit_dist)+"</div>":"")+
    '<div class="note">Stored coordinates are 0-based half-open (UCSC/BED). Displayed as 1-based inclusive.</div>';
}
function aliasTable(al){
  const by={}; al.forEach(a=>{(by[a.alias_type]=by[a.alias_type]||[]).push(a);});
  return "<table><tr><th>type</th><th>aliases</th></tr>"+Object.keys(by).sort().map(t=>{
    const seen=new Set(),out=[];
    by[t].forEach(a=>{const k=a.alias+"|"+(a.assignment||"");if(seen.has(k))return;seen.add(k);
      out.push('<span class="mono">'+esc(a.alias)+"</span>"+
        (a.assignment?' <span class="badge">'+esc(a.assignment)+"</span>":"")+
        (a.is_current?"":' <span class="badge retired">retired</span>'));});
    return "<tr><td>"+esc(t)+"</td><td>"+out.join("<br>")+"</td></tr>";}).join("")+"</table>";
}

// ---- locus graphic ----
async function drawLocus(d,co,asm){
  // called unawaited from render(); a throw here would otherwise be an invisible
  // rejected promise, leaving the "rendering…" placeholder on screen forever.
  try{ await drawLocus_(d,co,asm); }
  catch(e){ console.error("drawLocus",e);
    const g=$("gfx"); if(g) g.innerHTML='<div class="note" style="color:#a33">graphic failed: '+
      esc(e&&e.message||String(e))+'</div>'; }
}
// Repeat classes drawn by default. Everything else (Simple_repeat, Low_complexity,
// Satellite, tRNA/rRNA/snRNA/srpRNA, Unknown) is IN the shard and drawn only when the
// "all repeat classes" box is ticked -- see REP_LOW_INFO and the #repall handler.
const REP_INTERSPERSED=new Set(["LINE","SINE","LTR","LTR?","DNA","DNA?","Retroposon","RC","SINE?","LINE?"]);
const REP_COL={LINE:"#b5651d",SINE:"#2e8b8b","SINE?":"#2e8b8b",LTR:"#3b6ea5","LTR?":"#3b6ea5",
  DNA:"#6b7f3a","DNA?":"#6b7f3a",Retroposon:"#8f4c7a",RC:"#8f4c7a"};
const repCol=c=>REP_COL[c]||"#999";

async function drawLocus_(d,co,asm){
  asm=asm||"hg38";
  const isHg=asm==="hg38";
  // L is the lane-label gutter. At L=62 only ~9 monospace chars fit and 81% of
  // gene labels were clipped (e.g. "LOC124905662" rendered as "OC124905662").
  // L=110 holds ~18; longer labels are ellipsised by laneLabel() with the full
  // string in a <title>, so nothing is silently truncated.
  const PAD=1000, W=1080, L=110, R=14;
  const w0=Math.max(0,co.start-PAD), w1=co.end+PAD, span=w1-w0;
  const x=p=>L+(Math.min(Math.max(p,w0),w1)-w0)/span*(W-L-R);
  // gEVE ORFs and HERVarium domains are stored in hg38 coordinates ONLY, so they
  // are omitted on t2t rather than drawn at wrong positions. Gene models are NOT:
  // the bundle is assembly-keyed, hg38 from GENCODE and t2t from hs1 RefSeq, each
  // in its own assembly's coordinates. Index by assembly, never assume hg38.
  const GMB=await geneModels(), gm=(GMB&&GMB[asm])||{}, gkey=d.uid;
  const tx=(gm[gkey]||[]).filter(t=>t.txEnd>w0&&t.txStart<w1);
  const lanes=[];
  const showAll=(typeof document!=="undefined"&&$("repall")&&$("repall").checked);
  const reps=(d.repeats||[]).filter(r=>r.assembly===asm&&r.end>w0&&r.start<w1
                                       &&(showAll||REP_INTERSPERSED.has(r.rep_class)));
  // locus_segment is Telescope-derived and hg38-only. On t2t the LTR-class rows from
  // locus_repeat stand in, which is why the lane label names its source.
  const segs=isHg?(d.segments||[]).filter(s=>s.end>w0&&s.start<w1):[];
  lanes.push({label:isHg?"segments":"locus extent",h:16,draw:()=>segs.length?segs.map(s=>{
      const c=s.segment_class==="ltr"?"var(--ltr)":"var(--int)";
      const a=x(s.start),b=x(s.end),txt=(s.repName||"")+(s.segment_class==="ltr"?" (LTR)":"");
      return rect(a,0,b-a,13,c)+
        lbl((a+b)/2,9.5,txt,"#fff",8,"middle",fits(b-a,txt,8));}).join("")
    :rect(x(co.start),0,x(co.end)-x(co.start),13,"#c8c8d0")+
     lbl((x(co.start)+x(co.end))/2,9.5,
         isHg?"locus extent (no segments stored)":"locus extent (t2t; segments are hg38-only)",
         "#555",8,"middle",true)});
  // RepeatMasker lane: both assemblies. Labels are length-aware and separated, same
  // rule as the domain lane -- a 300 bp Alu is ~4 px wide at this scale.
  if(reps.length)lanes.push({label:"RepeatMasker",h:15,draw:()=>{
      let last=-1e9;
      return reps.slice().sort((p,q)=>p.start-q.start).map(r=>{
        const a=x(r.start),b=x(r.end),txt=r.rep_name||"";
        const ok=fits(b-a,txt,7)&&a-last>3; if(ok)last=b;
        return rect(a,0,Math.max(1.5,b-a),12,repCol(r.rep_class))+
               lbl((a+b)/2,8.7,txt,"#fff",7,"middle",ok);
      }).join("");}});
  const orfs=isHg?(d.geve||[]).filter(o=>o.hg38_orf_end>w0&&o.hg38_orf_start<w1):[];
  if(orfs.length)lanes.push({label:"gEVE ORFs",h:15,draw:()=>orfs.map(o=>{
      const a=x(o.hg38_orf_start),b=x(o.hg38_orf_end),txt=o.orf_class||"ORF";
      return arrow(a,b,12,"var(--orf)",o.orf_strand)+
        lbl((a+b)/2,9,txt,"#fff",8,"middle",fits(b-a,txt,8));}).join("")});
  const doms=isHg?(d.domains||[]).filter(o=>o.hg38_end>w0&&o.hg38_start<w1):[];
  if(doms.length)lanes.push({label:"HERVarium domains",h:15,draw:()=>{
      // labels only where the box holds them AND no drawn label is within 3px
      let last=-1e9;
      return doms.slice().sort((p,q)=>p.hg38_start-q.hg38_start).map(o=>{
        const a=x(o.hg38_start),b=x(o.hg38_end),txt=o.gene||o.domain_desc||"";
        const ok=fits(b-a,txt,8)&&a-last>3; if(ok)last=b;
        return rect(a,0,Math.max(2,b-a),12,"var(--dom)")+lbl((a+b)/2,8.7,txt,"#fff",8,"middle",ok);
      }).join("");}});
  // collapse isoforms to one lane per gene: union of exons, widest tx extent.
  // Drawing 6 lanes of the same gene wastes vertical space and hides whether any exon is in view.
  const byGene=new Map();
  tx.forEach(t=>{const k=t.name2||t.name;
    if(!byGene.has(k))byGene.set(k,{k,n:0,a:t.txStart,b:t.txEnd,strand:t.strand,ex:[]});
    const g=byGene.get(k); g.n++; g.a=Math.min(g.a,t.txStart); g.b=Math.max(g.b,t.txEnd);
    const es=t.exonStarts||[],ee=t.exonEnds||[];
    for(let i=0;i<es.length;i++) if(ee[i]>w0&&es[i]<w1) g.ex.push([es[i],ee[i]]);});
  [...byGene.values()].sort((p,q)=>p.a-q.a).forEach(g=>{
    g.ex.sort((p,q)=>p[0]-q[0]);
    const mg=[]; for(const e of g.ex){const l=mg[mg.length-1];
      if(l&&e[0]<=l[1])l[1]=Math.max(l[1],e[1]);else mg.push([e[0],e[1]]);}
    lanes.push({label:g.k+(g.n>1?" ("+g.n+" tx)":""),h:15,draw:()=>{
      let s=line(x(g.a),6,x(g.b),6,"var(--gene)",1);
      for(const [a,b] of mg) s+=rect(x(a),1.5,Math.max(1.5,x(b)-x(a)),9,"var(--gene)");
      s+=lbl(x(Math.min(g.b,w1))+4,9.5,g.strand,"#666",9,"start",true);
      // y=9.5 put this ON the gene line at y=6, which struck through the text.
      // 13 clears the line and still sits inside the 15px lane.
      if(!mg.length) s+=lbl((x(Math.max(g.a,w0))+x(Math.min(g.b,w1)))/2,13,
          "intron only \u2014 no exon in window","#7a5c8f",7.5,"middle",true);
      return s;}});});
  let y=0,body="";
  lanes.forEach(ln=>{
    body+='<g transform="translate(0,'+y+')">'+
      laneLabel(L-6,10,ln.label,"#555",9.5,L-6)+ln.draw()+"</g>";
    y+=ln.h+5;});
  // locus extent guides + axis
  const guides=line(x(co.start),0,x(co.start),y,"#b9c6d4",1,"2,2")+
               line(x(co.end),0,x(co.end),y,"#b9c6d4",1,"2,2");
  let axis=line(L,y+4,W-R,y+4,"#999",1);
  const ticks=5;
  // edge ticks anchor inward: a centred label at the last tick overflows the viewBox
  // and is clipped by the browser (the right-hand coordinate showed as "4,043,77").
  for(let i=0;i<=ticks;i++){const p=w0+(w1-w0)*i/ticks;
    const an=i===0?"start":(i===ticks?"end":"middle");
    axis+=line(x(p),y+4,x(p),y+8,"#999",1)+
      lbl(x(p),y+19,Math.round(p).toLocaleString(),"#666",9,an,true);}
  axis+=lbl(L,y+33,asm+" "+co.chrom+"  ·  window "+(w1-w0).toLocaleString()+" bp  ·  locus "+
        (co.end-co.start).toLocaleString()+" bp","#666",9.5,"start",true);
  $("gfx").innerHTML='<svg viewBox="0 0 '+W+" "+(y+42)+'" width="'+W+'">'+guides+body+axis+"</svg>"+
    (tx.length?"":'<div class="note">no '+(isHg?"GENCODE":"RefSeq")+
      ' transcript in this window'+
      (Object.keys(gm).length?"":" (no "+asm+" gene models in bundle)")+"</div>")+
    (d.segments.length?"":'<div class="note">This locus has no stored RepeatMasker segments — '+
      "only telescope-origin loci carry them. Bar shows the merged locus extent.</div>");
}
const rect=(x,y,w,h,f)=>'<rect x="'+x+'" y="'+y+'" width="'+Math.max(1,w)+'" height="'+h+'" fill="'+f+'" rx="1.5"/>';
const line=(x1,y1,x2,y2,c,w,dash)=>'<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+
  '" stroke="'+c+'" stroke-width="'+w+'"'+(dash?' stroke-dasharray="'+dash+'"':"")+"/>";
// Lane label, right-anchored in the L-px gutter. Ellipsises rather than letting the
// browser clip at the viewBox edge, and keeps the full name in a <title> tooltip.
const laneLabel=(x,y,t,c,s,gutter)=>{
  const per=0.60*s, max=Math.max(3,Math.floor((gutter-4)/per)), full=String(t||"");
  const cut=full.length>max?full.slice(0,max-1)+"\u2026":full;
  return '<text x="'+x+'" y="'+y+'" fill="'+c+'" font-size="'+s+
    '" text-anchor="end" font-family="ui-monospace,Menlo,monospace">'+
    (cut===full?"":'<title>'+esc(full)+'</title>')+esc(cut)+"</text>";
};
const lbl=(x,y,t,c,s,a,show)=>show?'<text x="'+x+'" y="'+y+'" fill="'+c+'" font-size="'+s+
  '" text-anchor="'+a+'" font-family="ui-monospace,Menlo,monospace">'+esc(t)+"</text>":"";
// monospace at font-size s is ~0.60*s per char; require the box to hold the string with padding
const fits=(w,t,s)=>w>=(String(t||"").length*0.60*s)+6;
function arrow(x1,x2,h,f,strand){
  const w=x2-x1,t=Math.min(7,Math.max(2,w*0.25));
  if(w<5)return rect(x1,0,w,h,f);
  return strand==="-"
    ?'<path d="M'+(x1+t)+' 0 H'+x2+' V'+h+' H'+(x1+t)+' L'+x1+' '+h/2+' Z" fill="'+f+'"/>'
    :'<path d="M'+x1+' 0 H'+(x2-t)+' L'+x2+' '+h/2+' L'+(x2-t)+' '+h+' H'+x1+' Z" fill="'+f+'"/>';
}
