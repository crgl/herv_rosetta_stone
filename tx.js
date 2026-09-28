/* Transcription-evidence rendering: FANTOM5 TSS / splice-site status, Snaptron
   junction arcs, and the per-unit (TU) detail view.

   Loaded after detail.js and reuses its primitives (rect, line, lbl, arrow,
   laneLabel, fits, esc, fmt, $). Kept in its own file because the locus detail
   renderer was already 320 lines and the TU view is a second, independent page.

   --- Presentation rules this file exists to hold -------------------------------

   ABSENCE IS THREE-VALUED. "not assessable" (no FANTOM coverage), "assessed,
   none found", and "evidence present" are three different cells, never two. A
   dash means the element could not be assessed; an explicit "none" means it was
   looked at and nothing was there. Never render the two the same way.

   COUNTS ARE PRE-CAP. Arc panels state "showing N of M" from arcs.shown and
   arcs.n_total. len(arcs.jx) is the shipped subset, not the evidence -- writing
   the shown count as if it were the total would understate support by orders of
   magnitude on repeat-dense loci (max in catalogue: 82,984 junctions in window).

   CROSS VS WITHIN IS THE WHOLE POINT, AND THE DIRECTION IS COUNTER-INTUITIVE. On a
   chimeric unit, arcs connecting two different member elements show a single
   spliced transcript running through both, which argues for KEEPING the merge;
   arcs confined within one member show no such link, which makes the unit a SPLIT
   candidate. See VERDICT below. They are drawn in different
   colours and counted separately, and a unit where neither class exists says so
   ("unanchored") rather than showing a blank panel. */

const TX_LEVEL_COL={stringent:"var(--tss3)",robust:"var(--tss2)",permissive:"var(--tss1)"};
const TX_LEVEL_RANK={permissive:1,robust:2,stringent:3};
const JX_CROSS="var(--jx-cross)", JX_WITHIN="var(--jx-within)", JX_ANY="var(--jx-any)", JX_UNANCH="var(--jx-anti)";
/* packed-reference classes on the LOCUS page: `span` = the element sits entirely
   within this intron, `edge` = one end inside the element, `both` reuses the
   within-element blue since it means the same thing at locus scale. */
const JX_SPAN="var(--jx-span)", JX_EDGE="var(--jx-edge)";

/* One evidence cell: a coloured chip at the strictest supported tier, a dash when
   the element was never assessable, the word "none" when it was assessed clean. */
function txChip(level,assessable,anti){
  if(!assessable) return '<span class="txna" title="no FANTOM CAT coverage over this '+
    'element — absence of evidence cannot be distinguished from absence of signal">'+
    "not assessable</span>";
  if(!level) return '<span class="txnone" title="assessed and no evidence found">none</span>';
  const c=TX_LEVEL_COL[level]||"var(--mut)";
  return '<span class="txchip" style="background:'+c+'">'+esc(level)+"</span>"+
    (anti?' <span class="txanti" title="evidence on the strand opposite the element">'+
      "antisense</span>":"");
}

/* Direction glyph: sense evidence points with the element, antisense against it.
   Drawn from the element's own strand, which the upstream layer already resolved --
   never re-derived from coordinates here. */
function txDir(rec,strand){
  if(!rec||!rec.assessable) return "—";
  const s=(rec.tss_level||rec.sj_level)?1:0, a=(rec.tss_anti||rec.sj_anti)?1:0;
  if(!s&&!a) return '<span class="txnone">none</span>';
  const fwd=strand==="-"?"◀":"▶", rev=strand==="-"?"▶":"◀";
  return (s?'<span class="txsense" title="sense — same strand as the element">'+fwd+
    " sense</span>":"")+(s&&a?" ":"")+
    (a?'<span class="txanti" title="antisense — opposite strand">'+rev+" antisense</span>":"");
}

/* Per-assembly TSS / splice table. Both assemblies are shown even though Snaptron
   arcs are hg38-only, because the FANTOM layer covers both and a missing t2t row
   is itself informative. */
function f5Rows(d){
  const st=(d.fantom5||{}).stats||{}, asms=["hg38","t2t"].filter(a=>st[a]);
  if(!asms.length) return "";
  return '<div class="f5row"><b>FANTOM5 CAGE peaks</b> (per locus; same primary data as FANTOM CAT, so not independent evidence)'+
    asms.map(a=>{const s=st[a];
      if(s.no_data) return '<br><span class="mono">'+a+"</span> \u2014 not assessable (contig carries no peaks)";
      const any=s.body||s.fkb||s.up||s.anti;
      return '<br><span class="mono">'+a+"</span> \u2014 "+(any?
        "body "+s.body+" sense \u00b7 first kb "+s.fkb+" \u00b7 1 kb upstream "+s.up+" \u00b7 antisense "+s.anti+
        (s.body&&s.max_tpm!=null?" \u00b7 body max "+s.max_tpm+" TPM":"")+(s.body&&s.breadth!=null?" in "+s.breadth+" librar"+(s.breadth===1?"y":"ies")+" \u2265 1 TPM":"")+
        (s.top?" \u00b7 top: "+esc(s.top):""):"assessed, none");}).join("")+"</div>";
}
/* Where the nearest sense TSS sits relative to the locus, from the signed distance
   to the 5' end (negative = upstream; bedtools closest -D a). */
function tssWhere(dist,span){
  if(dist==null) return "\u2014";
  if(dist<0) return (-dist).toLocaleString()+" bp upstream";
  if(span==null||dist<span) return dist.toLocaleString()+" bp into body";
  return (dist-span+1).toLocaleString()+" bp past 3\u2032 end";
}
function txPanel(d){
  const tx=d.tx||{}, asms=["hg38","t2t"].filter(a=>tx[a]);
  if(!asms.length){ const f=f5Rows(d);
    return f?'<div class="panel"><h2>Transcription initiation</h2>'+f+"</div>":""; }
  const strand=(coordOf(d,"hg38")||coordOf(d,"t2t")||{}).strand;
  const tierTag=(t,lvl)=>lvl?"":' <span class="note" style="font-style:normal">('+esc(t||"permissive")+")</span>";
  const rows=asms.map(a=>{
    const r=tx[a], co=coordOf(d,a)||{}, span=co.end!=null?co.end-co.start:null;
    const na=!r.assessable;
    const anti=r.tss_anti_n?r.tss_anti_n+(r.tss_anti_t&&r.tss_anti_t!==r.tss_t?' <span class="note" style="font-style:normal">('+esc(r.tss_anti_t)+")</span>":""):"0";
    return "<tr><td class=mono>"+esc(a)+"</td>"+
      "<td>"+txChip(r.tss_level,r.assessable,false)+"</td>"+
      "<td>"+txDir({assessable:r.assessable,tss_level:r.tss_level,tss_anti:r.tss_anti},strand)+"</td>"+
      "<td class=mono>"+(na?"\u2014":fmt(r.tss_body_n)+" / "+fmt(r.tss_up_n))+"</td>"+
      "<td class=mono>"+(na?"\u2014":anti)+"</td>"+
      "<td>"+(na?"\u2014":tssWhere(r.tss5p_dist,span)+tierTag(r.tss_t,r.tss_level))+"</td>"+
      "<td>"+txChip(r.sj_level,r.assessable,r.sj_anti)+"</td>"+
      "<td class=mono>"+(na?"\u2014":fmt(r.sj_donor_n)+" / "+fmt(r.sj_acceptor_n))+"</td></tr>";
  }).join("");
  return '<div class="panel"><h2>Transcription initiation &amp; splicing (FANTOM CAT)</h2>'+
    '<table class="t"><tr><th>assembly<th>TSS tier<th>direction<th>sense TSS body / 1 kb up'+
    "<th>antisense TSS<th>nearest sense TSS to 5\u2032 end<th>splice tier<th>sense donor / acceptor</tr>"+rows+"</table>"+
    '<div class="note">Tier is the <b>strictest</b> FANTOM CAT tier with sense evidence in the body or 1 kb upstream '+
    "(permissive \u2192 robust \u2192 stringent; a stricter tier implies the looser ones). <b>Every number in a row is "+
    "at that row\u2019s tier</b>, or at permissive when there is no sense evidence (marked in brackets). The nearest "+
    "sense TSS can lie outside the locus. \u201cnot assessable\u201d means no FANTOM CAT coverage over the element \u2014 "+
    "distinct from \u201cnone\u201d, which means assessed and silent.</div>"+f5Rows(d)+"</div>";
}

/* Snaptron summary for the locus page: counts only. The arc graphic lives in the
   locus map lane, and the unit-level classes live on the TU page. */
function jxPanel(d){
  const a=d.arcs||{};
  if(!a.n_total) return "";
  return '<div class="panel"><h2>Snaptron srav3h junctions</h2>'+
    '<dl class="kv">'+kv("junctions in ±1 kb window",a.n_total)+
    kv("shown in map",a.shown+(a.n_total>a.shown?" (top by sample support)":""))+
    kv("max support",Math.max.apply(null,a.jx.map(j=>j[2])).toLocaleString()+" samples")+
    "</dl>"+
    (a.n_total>a.shown?'<div class="note">'+a.shown+" of "+a.n_total.toLocaleString()+
      " junctions are drawn, ranked by the number of samples supporting each. A static "+
      "bundle cannot carry all 1.43M junctions at samples_count ≥ 10.</div>":"")+
    '<div class="note">hg38 only — Snaptron srav3h is not available for t2t.</div></div>';
}

/* Arc lane for the locus map. Junction arcs are drawn as quadratic curves above
   the axis; height encodes nothing (span is already the x-extent), colour encodes
   class, and opacity scales with log support so deep junctions read first. */
/* ---------------------------------------------- packed junction reference

   Decodes the bit-packed arcs the non-debug bundle ships (see pack_junctions.py
   for the authoritative layout):

     bits 63..41 (23)  donor offset, signed, relative to the locus 5' end
     bit  40     ( 1)  1 = junction strand differs from the locus strand
     bits 39..16 (24)  acceptor offset, signed, same origin
     bits 15..8  ( 8)  sample count,    10^(code/32)
     bits  7..0  ( 8)  coverage/sample, 10^(code/32)

   The words arrive as decimal strings in JSON, not numbers: a 64-bit value
   exceeds Number.MAX_SAFE_INTEGER (2^53), so parsing one as a double silently
   corrupts the donor field. BigInt is used for the shifts and the result is
   narrowed to Number only after masking, where every field fits in 24 bits.

   Output shape deliberately matches the debug `arcs.jx` rows -- absolute
   [start, end, sc, strand, canonical, class] -- so arcLane draws either source
   with no branch. */
var PJX_LOG = 32;

/* Round to 2 significant figures, for displaying log-quantised values without
   implying more precision than the 8-bit code carries. */
function sig2(v){
  if(!v || !isFinite(v)) return 0;
  const m = Math.pow(10, Math.floor(Math.log10(Math.abs(v))) - 1);
  return Math.round(v/m)*m;
}

function pjxDecode(pjx, co){
  if(!pjx || !pjx.w || !pjx.w.length) return {jx:[], n_total:0, shown:0};
  const five = co.strand === "+" ? co.start : co.end;
  const sgn  = co.strand === "+" ? 1 : -1;
  // Layout (see pack_junctions.py): 22 donor | 1 strand | 22 acceptor |
  // 1 anchor-mappable | 2 spare | 8 sc | 8 cov.  Both offsets are 22 bits, so
  // one mask and one sign-extension threshold serve both.
  const M22 = (1n<<22n)-1n;
  const jx = pjx.w.map(raw=>{
    const w = BigInt(raw);
    let d = (w>>42n) & M22;
    let a = (w>>19n) & M22;
    if(d >= 1n<<21n) d -= 1n<<22n;   // sign-extend 22-bit
    if(a >= 1n<<21n) a -= 1n<<22n;   // sign-extend 22-bit
    const anti = Number((w>>41n) & 1n);
    // 1 = at least one outer 50bp anchor window is wholly inside a pm151 easy
    // region, i.e. a split read can actually be placed at this junction.
    const anch = Number((w>>18n) & 1n);
    const sc  = Math.round(Math.pow(10, Number((w>>8n) & 0xFFn)/PJX_LOG));
    const cov = Math.pow(10, Number(w & 0xFFn)/PJX_LOG);
    // undo the 5'-relative, strand-signed encoding
    let p = five + Number(d)*sgn, q = five + Number(a)*sgn;
    const lo = Math.min(p,q), hi = Math.max(p,q);
    const cls = (lo < co.start && hi > co.end) ? "span"
              : (lo >= co.start && hi <= co.end) ? "both" : "edge";
    const st = anti ? (co.strand === "+" ? "-" : "+") : co.strand;
    return [lo, hi, sc, st, 1, cls, cov, anti, anch];
  });
  return {jx:jx, n_total: pjx.n||jx.length, shown: jx.length, packed:true};
}

/* Absolute support scale, identical on every locus page: width and opacity both
   follow log10(samples) from the Snaptron floor (10): 1.0 px, 1k 2.0 px, 100k 3.0 px. */
function jxW(sc){ return 1.0+0.5*Math.max(0,Math.min(4.6,Math.log10(Math.max(sc,10))-1)); }
function jxO(sc){ return 0.42+0.5*Math.max(0,Math.min(1,(Math.log10(Math.max(sc,10))-1)/4)); }
/* Arc lane. Baseline at y=H. Sense arcs (junction on the locus strand) rise above
   it, antisense arcs hang below it, to a depth of HB. Colour = class, width and
   opacity = absolute log support, dashed = neither 50 bp anchor mappable (pm151).
   Each arc is a <g class="jx"> holding a wide transparent hit path and the drawn
   path, so thin arcs are easy to hover and click. */
function arcLane(arcs,x,w0,w1,H,lstrand,HB){
  const jx=(arcs&&arcs.jx||[]).filter(j=>j[1]>w0&&j[0]<w1);
  if(!jx.length) return "";
  HB=HB||0;
  return jx.slice().sort((p,q)=>p[2]-q[2]).map(j=>{
    const a=x(j[0]), b=x(j[1]), cls=j[5]||"any";
    const col=cls==="cross"?JX_CROSS:cls==="within"?JX_WITHIN:
              cls==="unanchored"?JX_UNANCH:
              cls==="both"?JX_WITHIN:cls==="span"?JX_SPAN:
              cls==="edge"?JX_EDGE:JX_ANY;
    const anti=j.length>7?!!j[7]:(lstrand==="+"||lstrand==="-")&&j[3]!==lstrand;
    const measured=j.length>8, unanch=measured&&!j[8];
    const below=anti&&HB>0;
    const reach=below?HB-3:H-3;
    const h=Math.min(reach,5+reach*0.75*Math.min(1,(b-a)/420));
    const cy=below?H+h*2:H-h*2;
    const d='M'+a+" "+H+" Q"+((a+b)/2)+" "+cy+" "+b+" "+H;
    const sc=j[2], wdt=jxW(sc), op=jxO(sc)*(cls==="span"?0.8:1);
    const tip="<b>"+esc(cls==="span"?"spanning":cls)+" junction</b><br>"+
      '<span class="k">intron:</span> '+j[0].toLocaleString()+"\u2013"+j[1].toLocaleString()+
      " ("+(j[1]-j[0]).toLocaleString()+" bp)<br>"+
      '<span class="k">samples:</span> '+(arcs&&arcs.packed?"~"+sig2(sc).toLocaleString():sc.toLocaleString())+
      (arcs&&arcs.packed?" (log-quantised, \u00b14%)":"")+"<br>"+
      '<span class="k">strand:</span> '+esc(j[3])+(anti?" \u2014 antisense to the locus (drawn below)":" \u2014 sense (drawn above)")+
      (j.length>6&&j[6]?'<br><span class="k">coverage:</span> '+(arcs&&arcs.packed?"~":"")+j[6].toFixed(1)+" reads/sample":"")+
      (measured?'<br><span class="k">anchors:</span> '+(unanch?"NEITHER 50 bp anchor mappable (pm151) \u2014 dashed":"\u22651 50 bp anchor mappable (pm151)"):"");
    return '<g class="jx" data-tip="'+esc(tip)+'"><path class="hit" d="'+d+'" fill="none" stroke="transparent" stroke-width="9"/>'+
      '<path class="vis" d="'+d+'" fill="none" stroke="'+col+'" stroke-width="'+wdt.toFixed(2)+'"'+
      (unanch?' stroke-dasharray="5,3"':"")+' stroke-linecap="round" opacity="'+op.toFixed(2)+'"/></g>';
  }).join("");
}

/* ------------------------------------------------------------------ TU view */

/* Verdict semantics -- the direction here is the whole point, and it is easy to
   invert by intuition, so it is spelled out:

     cross-element junctions  ->  KEEP MERGED. A junction with one end in member A
     and the other in member B means a single spliced transcript runs through both.
     They are co-transcribed, so the merged unit is the right object.

     within-element only      ->  SPLIT CANDIDATE. Each member splices only inside
     itself; nothing observed ties them into one transcript, so the merge rests on
     proximity alone.

   The upstream layer assigns the verdict this way (cross_element_jx_n > 0 ->
   keep_merged_jx_supported) and the tri-panel figure labels it the same. An
   earlier draft of this file had the two descriptions swapped; the render sweep
   caught it by asserting the verdict text against the counts, which is why that
   assertion stays in simtu.mjs. */
const VERDICT={
  split_candidate_jx_supported:{t:"split candidate — junction-supported",c:"var(--bad)",
    d:"Junctions are confined within individual member elements — nothing splices "+
      "from one member into another. No observed transcript ties the members "+
      "together, so the merge rests on proximity rather than shared transcription."},
  keep_merged_jx_supported:{t:"keep merged — junction-supported",c:"var(--ok)",
    d:"Cross-element junctions connect member elements directly: a spliced "+
      "transcript runs from one member into another, so they are co-transcribed and "+
      "the merged unit is the right object to keep."},
  jx_uninformative:{t:"junction-uninformative",c:"var(--neutral)",
    d:"Too little junction evidence over this unit to speak either way. Absence here "+
      "is absence of data, not evidence for merging."},
  not_chimeric:{t:"not chimeric",c:"var(--jx-any)",
    d:"Single-group unit — the split question does not arise."}
};

function verdictPanel(s,arcs){
  const v=VERDICT[s.split_verdict]||{t:s.split_verdict||"—",c:"var(--mut)",d:""};
  return '<div class="panel"><h2>Split verdict</h2>'+
    '<div class="vbox" style="border-left:4px solid '+v.c+'"><b style="color:'+v.c+'">'+
    esc(v.t)+"</b><div class=\"note\">"+esc(v.d)+"</div></div>"+
    '<dl class="kv">'+
    kv("member loci",s.n_mem_loci)+kv("member groups",s.n_mem_groups)+
    kv("cross-element junctions",s.cross_element_jx_n)+
    kv("max cross support",s.cross_max_sc==null?null:s.cross_max_sc+" samples")+
    kv("within-element junctions",s.within_element_jx_n)+
    kv("max within support",s.within_max_sc==null?null:s.within_max_sc+" samples")+
    "</dl>"+
    '<div class="note">A <b>cross-element</b> junction argues for keeping the merge '+
    "(one transcript spans two members); <b>within-element</b> junctions argue for "+
    "splitting (each member splices only inside itself). "+
    "Membership classes are the audited ones: <b>internal</b> "+
    "junctions have both ends inside the unit, <b>boundary</b> junctions have one end "+
    "inside, and <b>intronic-context</b> junctions are host-gene introns that "+
    "<i>contain</i> the unit (median span 36 kb against a 2.4 kb median unit) — "+
    "containment, not readthrough of the element.</div>"+
    /* Summary counts above are the upstream Snaptron layer's; the arcs on the map
       are recomputed from the junction table for this view. They usually agree.
       Where they do not, both numbers are shown -- silently displaying one while
       drawing the other is how a page ends up contradicting itself. */
    (arcs&&arcs.recount
      ? '<div class="note warn">Arc classification on the map is recomputed here and '+
        "gives <b>"+arcs.n_cross+" cross</b> / <b>"+arcs.n_within+" within</b>, "+
        "against "+arcs.recount[0]+" / "+arcs.recount[1]+" in the upstream summary "+
        "above. The upstream window padding and support filter are not fully "+
        "reconstructable, so treat the counts as approximately concordant and the "+
        "verdict as resting on the upstream figure.</div>"
      : "")+"</div>";
}

/* Member-element track with junction arcs drawn over it. This is the view the
   split decision is actually made from: every member element in one frame, with
   the arcs that either cross between them or stay inside one. */
function drawTU(t){
  const a=t.arcs||{}, mem=t.members||[], ext=t.extent;
  if(!mem.length&&!ext&&!a.n_total){ $("gfx").innerHTML='<div class="note">no stored '+
    "coordinates for this unit in either assembly — nothing to draw</div>"; return; }
  const W=1080, L=110, R=14, ARCH=64;
  const xs=mem.length?mem.map(m=>m[2]):(ext?[ext[0]]:[0]);
  const xe=mem.length?mem.map(m=>m[3]):(ext?[ext[1]]:[0]);
  const w0=a.win?a.win[0]:Math.min.apply(null,xs)-1000;
  const w1=a.win?a.win[1]:Math.max.apply(null,xe)+1000;
  const span=Math.max(1,w1-w0);
  const x=p=>L+(Math.min(Math.max(p,w0),w1)-w0)/span*(W-L-R);

  // group colours: distinct hues per member group so a chimeric unit reads at a glance
  const groups=[...new Set(mem.map(m=>m[5]).filter(Boolean))];
  const GC=["var(--cat1)","var(--cat2)","var(--cat3)","var(--cat4)","var(--cat5)","var(--cat6)"];
  const gcol=g=>GC[Math.max(0,groups.indexOf(g))%GC.length];

  let body="", y=0;
  // arc lane on top, so arcs visually sit above the elements they connect
  if(a.jx&&a.jx.length){
    body+='<g transform="translate(0,'+y+')">'+laneLabel(L-6,ARCH,"junctions","var(--mut)",9.5,L-6)+
      arcLane(a,x,w0,w1,ARCH)+"</g>";
    y+=ARCH+6;
  }
  // No member mapping exists upstream: draw the unit extent as one unmapped bar.
  // Grey and explicitly labelled, so it cannot be mistaken for a resolved element.
  if(!mem.length&&ext){
    body+='<g transform="translate(0,'+y+')">'+
      laneLabel(L-6,11,"unit extent","var(--mut)",9.5,L-6)+
      rect(x(ext[0]),0,Math.max(1,x(ext[1])-x(ext[0])),14,"var(--rep-low)")+
      lbl((x(ext[0])+x(ext[1]))/2,10,"no member mapping","var(--ink)",8,"middle",
          fits(x(ext[1])-x(ext[0]),"no member mapping",8))+"</g>";
    y+=19;
  }
  // member elements, one lane each, labelled with group and combined_id
  mem.slice().sort((p,q)=>p[2]-q[2]).forEach(m=>{
    const [uid,cid,st,en,strand,grp]=m;
    const A=x(st), B=x(en), txt=(grp||"")+" "+(strand||"");
    body+='<g transform="translate(0,'+y+')">'+
      laneLabel(L-6,11,grp||uid,"var(--mut)",9.5,L-6)+
      arrow(A,0,B,14,gcol(grp),strand).replace("M"+A,"M"+A)+
      lbl((A+B)/2,10,txt,"var(--on-strong)",8,"middle",fits(B-A,txt,8))+
      '<title>'+esc(cid||uid)+"  "+esc(grp||"")+"  "+st.toLocaleString()+"–"+
      en.toLocaleString()+" ("+esc(strand)+")</title></g>";
    y+=19;
  });
  // axis
  let axis=line(L,y+4,W-R,y+4,"var(--axis)",1);
  for(let i=0;i<=5;i++){const p=w0+span*i/5;
    const an=i===0?"start":(i===5?"end":"middle");
    axis+=line(x(p),y+4,x(p),y+8,"var(--axis)",1)+
      lbl(x(p),y+19,Math.round(p).toLocaleString(),"var(--mut)",9,an,true);}
  axis+=lbl(L,y+33,esc(t.map_asm||"hg38")+" "+esc(t.map_chrom||a.chrom||"")+
    "  ·  window "+span.toLocaleString()+
    " bp  ·  "+mem.length+" member element"+(mem.length===1?"":"s"),"var(--mut)",9.5,"start",true);

  const leg='<div class="note"><span style="color:'+JX_CROSS+'">━</span> cross-element  '+
    '<span style="color:'+JX_WITHIN+'">━</span> within-element  '+
    '<span style="color:'+JX_UNANCH+'">━</span> unanchored (neither end inside a member)  '+
    "· arc opacity scales with log sample support</div>";
  // Two different reasons the shown set is smaller than the window total, and
  // they must not be described the same way. Cap-limited: the budget bit, and the
  // omitted arcs are lower-support. Anchoring-limited: most window junctions have
  // an end outside every member element, so they were never candidates -- calling
  // that "top 25 by support" would misdescribe the selection.
  let cap="";
  if(a.n_total>a.shown){
    const anch=a.n_anchored, lim=(anch!=null&&anch>a.shown);
    cap='<div class="note">Showing '+a.shown+" of "+a.n_total.toLocaleString()+
      " junctions in window. "+(a.n_cross!=null
        ? (lim?"Up to 25 cross-element and 25 within-element by sample support, so "+
             "support ranking cannot bury the split-relevant class."
             :"Only "+anch+" of "+a.n_total.toLocaleString()+" have both ends inside a "+
              "member element; the rest are not candidates for either class, so this "+
              "is all the anchored evidence, not a support cut-off.")
        : "Top by sample support.")+"</div>";
  }
  const unan=a.unanchored?'<div class="note">No junction in this window has '+
    "<b>both</b> ends inside a member element; the arcs shown start or land in the "+
    "flanks. The cross-element class is genuinely empty here — that is a real "+
    "observation, not missing data.</div>":"";
  $("gfx").innerHTML='<svg viewBox="0 0 '+W+" "+(y+42)+'" width="'+W+'">'+body+axis+
    "</svg>"+leg+cap+unan;
}

function renderTU(t,h){
  const s=t.snaptron||{}, tx=t.tx||{}, mem=t.members||[];
  const memRows=mem.length
    ? '<table class="t"><tr><th>group<th>combined_id<th>locus_uid<th>coordinates<th>strand</tr>'+
      mem.slice().sort((p,q)=>p[2]-q[2]).map(m=>
        "<tr><td>"+esc(m[5]||"—")+"</td><td class=mono>"+esc(m[1]||"")+
        '</td><td class=mono><a href="#" data-tuu="'+esc(m[0])+'">'+esc(m[0])+
        "</a></td><td class=mono>"+m[2].toLocaleString()+"–"+m[3].toLocaleString()+
        "</td><td class=mono>"+esc(m[4])+"</td></tr>").join("")+"</table>"
    : '<div class="note">no member coordinates in hg38</div>';
  const asms=["hg38","t2t"].filter(a=>tx[a]);
  const txRows=asms.map(a=>{const r=tx[a];
    return "<tr><td class=mono>"+esc(a)+"</td><td>"+txChip(r.tss_level,r.assessable,r.tss_anti)+
      "</td><td class=mono>"+fmt(r.tss_body_n)+" / "+fmt(r.tss_up_n)+"</td><td>"+
      txChip(r.sj_level,r.assessable,r.sj_anti)+"</td><td class=mono>"+fmt(r.sj_donor_n)+
      " / "+fmt(r.sj_acceptor_n)+"</td></tr>";}).join("");

  $("view").innerHTML=
    '<div class="panel"><div class="idline"><span class="cid">'+esc(t.tu_id)+"</span>"+
    '<span class="badge">transcriptional unit v0.2</span>'+
    (s.is_chimeric?'<span class="badge ambig">chimeric · '+fmt(s.n_mem_groups)+
      " groups</span>":"")+"</div>"+
    '<div class="note">Transcriptional units are the merge layer over catalogue loci; '+
    "this page is the internal view used to adjudicate whether a unit should be split. "+
    "Junction evidence is Snaptron srav3h (hg38).</div></div>"+
    '<div class="panel"><h2>Unit map — members and junctions</h2><div id="gfx"></div></div>'+
    verdictPanel(s,t.arcs)+
    '<div class="panel"><h2>Member elements — '+mem.length+"</h2>"+memRows+"</div>"+
    (txRows?'<div class="panel"><h2>Transcription (FANTOM5, unit extent)</h2>'+
      '<table class="t"><tr><th>assembly<th>TSS<th>TSS n body/up<th>splice'+
      "<th>donor/acceptor n</tr>"+txRows+"</table></div>":"")+
    '<div class="panel"><h2>Junction classes</h2><dl class="kv">'+
      kv("internal (both ends inside)",s.jx_internal_n)+
      kv("max internal support",s.jx_internal_max_sc)+
      kv("boundary (one end inside)",s.jx_boundary_n)+
      kv("intronic context (contains unit)",s.jx_intronic_ctx_n)+
      kv("median containing-intron span",s.intronic_ctx_median_span==null?null:
         Number(s.intronic_ctx_median_span).toLocaleString()+" bp")+
      "</dl></div>";
  drawTU(t);
  // member links jump to the locus page, so a unit is one click from any element
  [...document.querySelectorAll("a[data-tuu]")].forEach(el=>el.onclick=async e=>{
    e.preventDefault();
    try{ const d=await getLocus(el.dataset.tuu);
      if(!d) throw new Error("locus "+el.dataset.tuu+" absent from shard");
      render(d,{k:el.dataset.tuu,kind:"u",type:"uid",uids:[el.dataset.tuu],cur:1});
    }catch(err){ showError("Loading locus",err); }});
}

/* TSS tick marks for the locus/TU map: a caret at the 5′ end when TSS evidence
   exists, coloured by tier and pointing in the direction of transcription. */
/* FANTOM CAT marks on the junction lane. The sense arrow is drawn AT the nearest
   sense TSS of the row's tier when that TSS falls in the drawn window; otherwise a
   small flag at the 5' end says evidence exists but its position is off-window. */
function tssMarks(rec,co,x,strand){
  if(!rec||!rec.assessable||!rec.tss_level||!(strand==="+"||strand==="-")) return "";
  const c=TX_LEVEL_COL[rec.tss_level]||"var(--mut)", dd=strand==="-"?-1:1;
  const span=co.end-co.start, dist=rec.tss5p_dist;
  const pos=dist==null?null:(strand==="+"?co.start+dist:co.end-1-dist);
  const onmap=pos!=null&&pos>=co.start-1000&&pos<co.end+1000;
  const where=tssWhere(dist,span);
  const X=x(onmap?pos:(strand==="-"?co.end:co.start));
  const tip="<b>FANTOM CAT sense TSS</b><br>"+'<span class="k">tier:</span> '+esc(rec.tss_level)+
    (rec.tss_body_n+rec.tss_up_n?"<br>"+'<span class="k">in body / 1 kb up:</span> '+rec.tss_body_n+" / "+rec.tss_up_n:"")+
    "<br>"+'<span class="k">nearest to 5\u2032 end:</span> '+where+
    (onmap?"":"<br>(off the map window \u2014 this mark is a flag at the 5\u2032 end, not the TSS position)")+
    (rec.tss5p_score!=null?"<br>"+'<span class="k">CAGE score:</span> '+rec.tss5p_score.toFixed(1):"");
  return '<path d="M'+X+' 12 L'+X+' 2 L'+(X+7*dd)+' 2 L'+(X+4*dd)+' 5" fill="none" stroke="'+c+'" stroke-width="1.8"'+
    (onmap?"":' stroke-dasharray="2,1.5"')+' data-tip="'+esc(tip)+'"></path>';
}
